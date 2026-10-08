"""Persist run state before submitting a cloud request."""
import time
from .common import write_json, fingerprint, log

class StateWriter:
    """Persist reservations locally before throttled, retryable GCS uploads.

    The coordinator's exclusive lock makes repeating the same state upload
    idempotent. Keep local progress even if all cloud attempts fail; resume
    reconciles reserved/accepted jobs with Batch before submitting anything.
    """
    def __init__(self, gcs, path, uri, *, clock=time.monotonic, sleep=time.sleep):
        self.gcs, self.path, self.uri = gcs, path, uri
        self.clock, self.sleep = clock, sleep
        self.next_upload = 0.
        self.saved = None

    def save(self, state):
        write_json(self.path, state)
        digest = fingerprint(state)
        if digest == self.saved:
            return
        for attempt in range(8):
            delay = self.next_upload-self.clock()
            if delay > 0:
                self.sleep(delay)
            try:
                self.gcs.state(self.path, self.uri)
            except Exception as error:
                # google.api_core HTTP errors expose their numeric status code.
                # Permanent/authentication errors must remain visible.
                if getattr(error, 'code', None) not in (408, 429, 500, 502, 503, 504) or attempt == 7:
                    raise
                log(f'State upload returned HTTP {error.code}; retrying saved progress')
                self.sleep(min(2**(attempt+1), 30))
            else:
                self.saved = digest
                return
            finally:
                # Cloud Storage permits one write/second to an object name.
                self.next_upload = self.clock()+1.1
