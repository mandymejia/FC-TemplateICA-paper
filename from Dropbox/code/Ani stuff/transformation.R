
##### The transformation function

  transform.s = function(s,rows){
  		if(!rows){s = t(s)}
 		sddev= apply(s,1,sd)
		for(k in 1:dim(s)[1]){s[k,]=(s[k,]-mean(s[k,]))/(sddev[k])}
		sthird=apply(s^3,1,mean)
		ssign=sign(sthird)
		for(k in 1:dim(s)[1]){s[k,]=s[k,]*ssign[k]}
		sthird=sthird*ssign
		or=order(sthird)
		s=s[or,]

		return(s)	
  }

  transform.sa = function(s,a,rows){
   		if(!rows){
   			s = t(s)
   			a=t(a)}
 		sddev= apply(s,1,sd)
		for(k in 1:dim(s)[1]){s[k,]=(s[k,]-mean(s[k,]))/(sddev[k])}
		sthird=apply(s^3,1,mean)
		ssign=sign(sthird)
		for(k in 1:dim(s)[1]){s[k,]=s[k,]*ssign[k]}
		sthird=sthird*ssign
		or=order(sthird)
		s=s[or,]
		for(k in 1:dim(s)[1]){a[,k]=a[,k]*(sddev[k])*ssign[k]}
		a=a[,or]
		return(c(list(s),list(a)))	
  }

