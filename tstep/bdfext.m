i=sqrt(-1.); th=0:10000; th=2*pi*th'/10000;
em1=exp(-1*i*th); em2=exp(-2*i*th); em3=exp(-3*i*th); e0=1. + 0*em1;
ep1=exp( 1*i*th); ep2=exp( 2*i*th); ep3=exp( 3*i*th);
hold off
xm=-1.4;xM=7.4;ym=-4.4;yM=4.4;xa=[xm xM];y0=0*xa;ya=[ym yM];x0=0*ya;
plot(xa,y0,'k-',x0,ya,'k-'); axis equal; axis([xm xM ym yM]); hold on;
bdf1=(1-1*em1)/1;
bdf2=(3-4*em1+em2)/2;
bdf3=(11-18*em1+9*em2-2*em3)/6;
text( 1.65,0.98,'k=1','FontSize',18);
text( 3.25,1.97,'2','FontSize',18); text( 4.63,3.13,'3','FontSize',18);
text(-1.09,4.03,'Im $$\lambda \Deltat$$','FontSize',18);
text( 5.60,0.40,'Re $$\lambda \Deltat$$','FontSize',18);
title('BDFk Neutral Stability Curve','FontSize',18);
hold off
xm=-2.4;xM=0.4;ym=-1.4;yM=1.4;xa=[xm xM];y0=0*xa;ya=[ym yM];x0=0*ya;
plot(xa,y0,'k-',x0,ya,'k-'); axis equal; axis([xm xM ym yM]); hold on;

bdfext1= bdf1./(em1);              plot(bdfext1,'g-');
bdfext2= bdf2./(2*em1-em2);        plot(bdfext2,'r-');
bdfext3= bdf3./(3*em1-3*em2+em3);  plot(bdfext3,'b-');

title('BDF/EXtk Neutral Stability Curve','FontSize',18);
text(-1.65,0.98,'k=1','FontSize',18);
text(-1.10,0.72,'2','FontSize',18); text(-0.83,0.53,'3','FontSize',18);
xlabel('Im $$abc\lambda \Deltat$$');
ylabel('Re $$abc\lambda \Deltat$$');
