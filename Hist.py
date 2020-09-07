
class Histogram:
      def __init__(self,min_lim,max_lim,width_bin,size):
          self.min = min_lim
          self.max = max_lim
          self.width_bin = width_bin
          self.size = size

      def idx_hist(self,value):
          return round((value - self.min)/self.width_bin)
