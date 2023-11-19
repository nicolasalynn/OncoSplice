from dataclasses import dataclass

@dataclass
class SpliceSite(object):
    pos: int
    ss_type: int
    prob: float

    def __post_init__(self):
        pass

    def __lt__(self, other):
        return self.pos < other.pos

