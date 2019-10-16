png('intervention_interaction_diagram.png', width=1600, height=1200, res=180)
par(mar=c(2.1,2.1,0,0))
x = seq(-3,30,0.1);
x2= seq(-3,31,0.1);
shape=7;
shape2=20;

y = function(v,s) { (atan(v)+pi/2)**s }

y2_rescale = y(max(x)-10,shape)/y(max(x),shape2)
y2 = function(v,s) { y(v+1,shape2)*y2_rescale }

toX = function(v) { v*(max(x) - min(x)) + min(x) }
#toY = function(v) { v*(y(max(v)) - y(min(v))) - y(min(v)) }

#plot(x,y(x,shape)(atan(x)+pi/2)**shape, ylim=c(0,1.05*pi**shape), type='l', axes=F, xlab='', ylab='', yaxs='i', xaxs='i', lwd=3)
plot(x,y(x,shape), ylim=c(0,1.05*pi**shape), type='l', axes=F, xlab='', ylab='', yaxs='i', xaxs='i', lwd=3, col='red')
#lines(x2-1, y(x2,shape2)*y2_rescale, type='l', col='blue', lwd=3)
lines(x2, y2(x2), type='l', col='blue', lwd=3)

title(xlab=expression(paste('Effective reproductive number ', italic(R))), line=1)
title(ylab='Number of cases per year', line=1);
abline(v=0.5, lty=3);
#abline(h=pi**shape, lty=3);
abline(h=y(70,shape), lty=3);
axis(1, labels=F, at=range(x));
axis(2, labels=F, at=c(0,pi**shape));

text(max(x)*0.8, y(50,shape), pos=3, labels='Maximum possible epidemic size')
#text(y(1,shape), pi**shape, pos=3, labels='Maximum possible epidemic size')
text(0, 0.2*(pi**shape), srt=90, labels='Epidemic threshold')

addArrow = function(.x, .yfn, .col, .shape, .y_rescale=1) {
    print(.x)
    lines(toX(.x[c(1,1)]), .yfn(toX(.x[1:2]), .shape)*.y_rescale, col=.col, lwd=1.5)

    arrows(
      toX(.x[1]),
    .yfn(toX(.x[2]), .shape)*.y_rescale,
      toX(.x[2]),
    .yfn(toX(.x[2]), .shape)*.y_rescale,
    length=0.2,
    lwd=1.5, col=.col)
}

init = 0.94
#lines(toX(c(init/2, init, init)), y(toX(c(init/2, init/2, init)),shape))
#lines(toX(c(init/4, init/2, init/2)), y(toX(c(init/4, init/4, init/2)),shape))

addArrow(c(init, init/2), y, 'red', shape)
addArrow(c(init/2, init/4), y, 'red', shape)

points(toX(c(init, init/2, init/4)), y(toX(c(init, init/2, init/4)), shape), pch=19, cex=1.5, col='red')
#text(toX(c(init, init/2, init/4)), 1.02*y(toX(c(init, init/2, init/4)), shape), pos=2, labels=c('A', 'B', 'C')) 
text(toX(c(init, init/2, init/4)), 1.02*y(toX(c(init, init/2, init/4)), shape), pos=2, labels=c('No intervention', 'With 1\nintervention', '(B) Interventions\nprevent\ninfection'), col='red') 


reduction = y(toX(init/2),shape)/y(toX(init),shape)
lines(toX(c(init, init)), reduction*y(toX(c(init, init/2)),shape), col='red', lty=2)

points(toX(c(init,init)), c(y(toX(init/2),shape), reduction*y(toX(init/2),shape)), cex=1.5, col='red', bg=c(NA, 'white'), pch=21)
text(toX(init), reduction*y(toX(init/2),shape), labels='\n\n(A) Interventions\nonly prevent\nsymptoms', pos=2, col='red')


init2 = 0.42
#lines(toX(c(init2/2, init2, init2))-1, y(toX(c(init2/2, init2/2, init2)),shape2)*y2_rescale, col='blue')
#lines(toX(c(init2, init2))-1, y(toX(c(init2, init2/2)),shape2)*y2_rescale, col='blue')
#arrows(
#  toX(c(init2))-1,
#y(toX(c(init2/2)),shape2)*y2_rescale,
#  toX(c(init2/2))-1,
#y(toX(c(init2/2)),shape2)*y2_rescale,
#lwd=1.5, col='blue')

addArrow(c(init2, init2/2), y2, 'blue', shape2)

addArrow(c(init2/2, init2/4), y2, 'blue', shape2)

#lines(toX(c(init2/4, init2/2, init2/2))-1, y(toX(c(init2/4, init2/4, init2/2)),shape2)*y2_rescale, col='blue')
points(toX(c(init2, init2/2, init2/4)), y2(toX(c(init2, init2/2, init2/4)), shape2), pch=19, col='blue', cex=1.5) 

#text(toX(init2)-1, 0.98*y(toX(init2), shape2)*y2_rescale, pos=4, labels='A', col='blue') 
#text(toX(init2/2)-1, 1.02*y(toX(init2/2), shape2)*y2_rescale, pos=2, labels='B', col='blue') 
text(toX(init2/4), y2(toX(0.01+init2/4), shape2), pos=2, labels='(C)', col='blue') 

text(toX(0.5), y(toX(0.12), shape), labels='Red/blue indicates disease system\n\nEach infection-preventing intervention (solid circles)\nreduce transmisison by 50%', pos=4)
dev.off()
