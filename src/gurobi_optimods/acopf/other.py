'''

      if  minanglerad >= -.5*math.pi:
	lbound = maxprod*sin(minanglerad)


else if  minanglerad >= - math.pi:
  lbound = -maxprod

else if  minanglerad >= - 1.5*math.pi:
  ubound = maxprod*max( sin(maxanglerad), sin(minanglerad))
	lbound = -maxprod
      }
      else{
	ubound = maxprod
	lbound = -maxprod 
      }
    }
    else if(maxanglerad <= math.pi){
      ubound = maxprod

      if (minanglerad >= -.5*math.pi){
	lbound = maxprod*sin(minanglerad)
      }
      else if (minanglerad >= - math.pi){
	lbound = -maxprod
      }
      else if (minanglerad >= - 1.5*math.pi){
	lbound = -maxprod
      }
      else{
	lbound = -maxprod 
      }
    }
    else if(maxanglerad <= 1.5*math.pi){
      ubound = maxprod

      if (minanglerad >= -.5*math.pi){
	lbound = maxprod*min(sin(maxanglerad), sin(minanglerad))
      }
      else if (minanglerad >= - math.pi){
	lbound = -maxprod
      }
      else if (minanglerad >= - 1.5*math.pi){
	lbound = -maxprod
      }
      else{
	lbound = -maxprod 
      }

    }
else:
  ubound = maxprod
      lbound = -maxprod

'''
