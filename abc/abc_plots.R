require("RSQLite")
require(beanplot)
drv = dbDriver("SQLite")
db = dbConnect(drv, "refit.sqlite")
abc = dbGetQuery(db, 'select J.*, P.*, M.* from jobs J, parameters P, metrics M where J.serial = P.serial and J.serial = M.serial')
extra_serials = which(names(abc)== 'serial')[-1] 
abc = abc[,-c(extra_serials)]
abc = subset(abc, select=-c(startTime, duration, attempts, seed))
abc$post_bool = abc$posterior >= 0

par_cols = 6:11
met_cols = c(12, 14:24)
incomplete_sets = unique(abc$smcSet[abc$status!='D']) # 6
all_sets = unique(abc$smcSet)
complete_sets = setdiff(all_sets, incomplete_sets)
last_complete_set = max(complete_sets)

pdf('marginal_pars.pdf', width=8.5, height=11)
par(mfrow=c(6,1))
par(mar=c(2.1, 4.1, 1.1, 0.5))
for (col in par_cols) {
    #prior = hist(abc[abc$smcSet==last_complete_set, col], nclass=20, plot=F, freq=F)
    #post  = hist(abc[abc$smcSet==last_complete_set & abc$posterior >= 0, col],nclass=20, plot=F, freq=F)
    colname = names(abc)[col];
    cat(colname)
    #boxplot( abc[,colname] ~ abc$post_bool*abc$smcSet, col=(c('black', 'red')), main=colname)
    #beanplot( abc[,colname] ~ abc$post_bool*abc$smcSet, col=(c('black', 'red')), main=colname)
    beanplot( abc[,colname] ~ abc$post_bool*abc$smcSet, what=c(1,1,1,0), col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1,1, grey(0.7))), main='', ylab=colname, side='both')
    #prior = density(abc[abc$smcSet==last_complete_set, col])
    #post  = density(abc[abc$smcSet==last_complete_set & abc$posterior >= 0, col])
    #plot(post, col='red', main='', xlab=names(abc)[col], lwd=2)
    #plot(post, col=rgb(1,0,0,1/4), add=T)
    #lines(prior, col='black', lwd=2)
}
dev.off()

ylims = list(c(0,800), c(0,5), c(0,200), c(0,350), c(0,600), c(10,30000), c(10e-1, 10e4), c(-2, 6), c(0, 1), c(0, 1), c(0, 1), c(-6, 4), c(-0.1, 0.1))

pdf('marginal_mets.pdf', width=8.5, height=22)
par(mfrow=c(12,1))
par(mar=c(2.1, 4.1, 1.1, 0.5))
for (col in met_cols) {
    colname = names(abc)[col];
    cat(paste(colname, '\n'))
    done = complete.cases(abc[,colname])
    if (colname %in% c('beta0', 'beta1')) {
        done = abc[,colname] < 100 & done 
    }
    d = abc[done, colname]
    #pl_ymin = max(min(d), mean(d) - 4*IQR(d))
    #pl_ymax = min(max(d), mean(d) + 4*IQR(d))
    beanplot( d ~ abc$post_bool[done]*abc$smcSet[done], what=c(1,1,1,0), 
              col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1,1, grey(0.7))), 
              main='', ylab=colname, side='both', ylim=unlist(ylims[col-max(par_cols)]) )
}

dev.off()

#pdf('marginal_pars.pdf', width=8.5, height=11)
#par(mfrow=c(6,1))
#par(mar=c(4.1, 4.1, 0.5, 0.5))
#
#tmp=abc[ , c(par_cols)];
##tmp$smcSet = abc$smcSet;
#tmp$post = abc$posterior >= 0
#
#boxplot( post ~. , data= tmp);
#dev.off()
