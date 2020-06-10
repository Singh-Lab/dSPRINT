import tempfile
from threading import Lock


class TwoStageRunner:
    def __init__(self, stage1_file=None, stage2_file=None):
        self.lock = Lock()
        self.stage1_file = stage1_file or next(tempfile._get_candidate_names())
        self.stage2_file = stage2_file or next(tempfile._get_candidate_names())

    def stage1(self, *args, **kwargs):
        with self.lock:
            self._stage1(*args, **kwargs)
            with open(self.stage1_file, 'w'):
                pass

    def _stage1(self, *args, **kwargs):
        pass

    def stage2(self, *args, **kwargs):
        with self.lock:
            self._stage2(*args, **kwargs)
            with open(self.stage2_file, 'w'):
                pass

    def _stage2(self, *args, **kwargs):
        pass
