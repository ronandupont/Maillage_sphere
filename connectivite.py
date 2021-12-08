from numpy import *
def f(nb_pts):
   a=log((nb_pts-2)/10)
   print(a)
   a/=log(4)
   print(a)
   return(a)

print(f(13))
