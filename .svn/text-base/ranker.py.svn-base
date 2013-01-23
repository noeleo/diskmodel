"""
This class stores the rankings of the lowest minimums.
"""
class Ranker:
    
    def __init__(self, numToTrack):
        self.best = {}
        self.numToTrack = numToTrack
        self.benchmark = float('inf')
    
    def update(self, item, score):
        if score <= self.benchmark:
            # add this one to the ranking
            self.best[score] = item
            if len(self.best) > self.numToTrack:
                # we need to replace one of the rankings
                self.best.pop(self.benchmark)
                self.benchmark = max(self.best.keys())
            elif len(self.best) == self.numToTrack:
                # this only happens once when we need to
                # to set a benchmark
                self.benchmark = max(self.best.keys())
    
    def getRankings(self):
        return self.best.copy()