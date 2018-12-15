
   if verbose:
         print 'Loaded random catalogue....t =',time()-self.s,'sec'
         print rantop-1, 'points in random catalogue'

         if flag == True:
            print ' '
            print 'WARNING: Fewer random points than data points.  This'
            print 'might lead to statistical imprecision in Landy-Szalay'
            print 'estimator.  INCREASE SIZE OF INPUT RANDOM CATALOGUE!!'
            print ' '
            while True:
               x = raw_input("Continue Anyway? (y/n)>: ")
               if x == 'n' or x == 'N':
                  return
               elif x == 'y' or x == 'yes':
                  break
               sleep(1)
