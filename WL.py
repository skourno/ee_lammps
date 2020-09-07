import numpy as np

class WL_histogram:
      def __init__(self,n_sub,lnf,ratio_acc,lnf_crit):
          self.lnf = lnf
          self.ratio_acc = ratio_acc
          self.lnf_crit = lnf_crit
          self.n_sub = n_sub
          self.visits = np.zeros(n_sub,np.int64)
          self.wts    = np.ones(n_sub,np.double)

      def reset_visits(self):
          self.visits[:] = 0

      def incr_visits(self,i_sub):
          self.visits[i_sub] += 1

      def penalize(self,i_sub):
          self.incr_visits(i_sub)
          self.wts[i_sub] -= self.lnf
          self.check_visits()

      def halve_penalty(self):
          self.lnf *= 0.5
          self.reset_visits()

      def check_visits(self):
          ratio_current = np.min(self.visits/np.mean(self.visits))
          if (ratio_current > self.ratio_acc):
             self.halve_penalty()

      def write_to_a_file(self,file_WL):
          file_WL.write('%6.5f \n' %(self.lnf))
          file_WL.write("\n")
          mean_visits = np.mean(self.visits)
          for i_sub in range(self.n_sub):
             file_WL.write('%6.5f %6.5f %6.5f \n' %(i_sub*0.1, self.visits[i_sub]/mean_visits, \
                       self.wts[i_sub]))
          file_WL.write("\n")
          file_WL.flush()
