"""Private-workspace GCS storage and authoritative Batch status queries."""
import json
from .common import sha256

class GCS:
    def __init__(self,project):
        from google.cloud import storage
        self.client=storage.Client(project=project)
        self.project=project

    def blob(self,uri):
        bucket,key=uri.removeprefix('gs://').split('/',1)
        return self.client.bucket(bucket,user_project=self.project).blob(key)

    def put(self,path,uri):
        blob=self.blob(uri)
        digest=sha256(path)
        if blob.exists():
            blob.reload()
            if blob.metadata and blob.metadata.get('sha256')==digest:
                return
            raise ValueError('Refusing to overwrite a different frozen artifact: '+uri)
        blob.metadata={'sha256':digest}
        blob.upload_from_filename(str(path),if_generation_match=0,timeout=3600)

    def json(self,uri):
        b=self.blob(uri)
        return json.loads(b.download_as_bytes()) if b.exists() else None

    def state(self,path,uri):
        # Mutable progress record, separate from immutable input/result artifacts.
        self.blob(uri).upload_from_filename(str(path))

    def metadata(self,uri):
        b=self.blob(uri);b.reload()
        return {'generation':b.generation,'size':b.size,'crc32c':b.crc32c}

    def objects(self,uri):
        bucket,prefix=uri.removeprefix('gs://').split('/',1)
        return self.client.list_blobs(self.client.bucket(bucket,user_project=self.project),prefix=prefix)

    def completed(self,uri,manifest_id):
        result=self.json(uri+'/COMPLETE.json')
        if result is None:return None
        if result['manifest_id']!=manifest_id:raise ValueError('Completion fingerprint mismatch: '+uri)
        root=uri.removeprefix('gs://').split('/',1)[1]+'/'
        actual={b.name[len(root):]:b.size for b in self.objects(uri+'/')}
        for name,record in result['files'].items():
            if actual.get(name)!=record['bytes']:
                raise ValueError('Incomplete delocalization: '+uri+'/'+name)
        return result


class StatusLookupError(RuntimeError):
    pass


def statuses(config,run_id):
    from google.cloud import batch_v1
    from google.api_core import exceptions
    from requests.exceptions import ConnectionError, Timeout
    client=batch_v1.BatchServiceClient(transport='rest')
    try:
        regions = config.get('batch_regions', [config['region']])
        if (not isinstance(regions, list) or not regions
                or any(not isinstance(region, str) or not region for region in regions)
                or len(set(regions)) != len(regions)):
            raise ValueError('Invalid Batch status regions')
        result = {}
        for region in regions:
            rows=client.list_jobs(request={
                'parent':f"projects/{config['project']}/locations/{region}",
                'filter':f'labels.lai-run="{run_id}"','page_size':1000},timeout=30)
            for j in rows:
                key=j.labels['job-id']
                if key in result:
                    raise ValueError('Duplicate Batch job identity across regions')
                result[key]={'name':j.name,'region':region,'labels':dict(j.labels),'status':{
                    'state':j.status.state.name,
                    'runDuration':str(j.status.run_duration.total_seconds())+'s',
                    'statusEvents':batch_v1.JobStatus.to_dict(j.status)['status_events']}}
        return result
    except (exceptions.ServiceUnavailable,exceptions.DeadlineExceeded,
            exceptions.TooManyRequests,exceptions.InternalServerError,
            ConnectionError,Timeout) as error:
        # The REST transport can expose requests errors directly, including
        # while fetching a later page. Discard the entire partial lookup and
        # let the coordinator retry without making scheduling decisions.
        raise StatusLookupError(str(error)) from error
    finally:
        client.transport.close()
