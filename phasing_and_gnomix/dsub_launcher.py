"""Workspace-compatible dsub 0.5.2 entry point, without editing installed dsub."""
from functools import wraps
from dsub.providers import google_utils, google_batch, google_batch_operations

google_utils.CLOUD_SDK_IMAGE = ('gcr.io/google.com/cloudsdktool/cloud-sdk@'
                              'sha256:c6c6bec5b8b94e4eb2f8ade1f17cf4cb01f07973b7fbb3a4843f0a61592ad5ed')
google_utils.PREPARE_CMD = google_utils.PREPARE_CMD.replace('| python -c','| python3 -c')
google_batch._CONTINUOUS_LOGGING_CMD = google_batch._CONTINUOUS_LOGGING_CMD.replace('| python -c','| python3 -c')


def compatible_instance_policy(builder):
    """Use a supported boot disk for M3 without changing attached data disks."""
    @wraps(builder)
    def build(*args, **kwargs):
        policy = builder(*args, **kwargs)
        if (policy.machine_type.startswith('m3-')
                and policy.boot_disk.type_ in ('', 'pd-standard')):
            policy.boot_disk.type_ = 'pd-balanced'
        return policy
    return build


google_batch_operations.build_instance_policy = compatible_instance_policy(
    google_batch_operations.build_instance_policy)

if __name__ == '__main__':
    from dsub.commands.dsub import main
    main()
