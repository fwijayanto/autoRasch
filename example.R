##### Example ######

raschdata.multidim.new <- response.dico.gen()

rasch.multidim.bckstep.new <- bck.stepwise.selection(X=raschdata.multidim.new,fix.set = c(1:18), tracked = TRUE, maxcore = 2, namefile = "temp.rasch.multidim.bckstep")


gpcm.multidim.bckstep.new <- bck.gpcm.stepwise.selection(X=pcmdata.multidim,fix.set = c(1:18), tracked = TRUE, maxcore = 2, namefile = "temp.gpcm.multidim.bckstep")

fshd.bckstep.new <- bck.gpcm.stepwise.selection(X=fshd.minus49patients.minus10items[,4:162],fix.set = c(1:149), tracked = TRUE, maxcore = 40, namefile = "temp.fshd.bckstep")

fshd.best.karlien.46.1pl <- jml.est.pcm(pspp.p.nona[,4:162][,c(94,73,135,45,97,75,61,77,87,124,86,47,123,38,67,83,72,70,106,110,66,91,112,63,7,1,118,116,114,29,53,119,102,98,111,48,71,4,59,137,23,99,26,8,101,20)], desc = Qmap[c(94,73,135,45,97,75,61,77,87,124,86,47,123,38,67,83,72,70,106,110,66,91,112,63,7,1,118,116,114,29,53,119,102,98,111,48,71,4,59,137,23,99,26,8,101,20),3])
fshd.best.karlien.46.fit <- itemfit.orm(fshd.best.karlien.46.1pl)

fshd.best.orm.46.1pl <- jml.est.pcm(pspp.p.nona[,4:162][,pspp.gpcm.backstep.el.1$item.1[46,1:46]], desc = Qmap[pspp.gpcm.backstep.el.1$item.1[46,1:46],3])
fshd.best.orm.46.fit <- itemfit.orm(fshd.best.orm.46.1pl)

fshd.best.karlien.1pl.ld <- cor(fshd.best.karlien.46.fit$trace.mat$resid.mat, method = "pearson", use = "pairwise.complete.obs")
fshd.best.karlien.ld <- mean(abs(fshd.best.karlien.1pl.ld[lower.tri(fshd.best.karlien.1pl.ld)]))

fshd.best.orm.1pl.ld <- cor(fshd.best.orm.46.fit$trace.mat$resid.mat, method = "pearson", use = "pairwise.complete.obs")
fshd.best.orm.ld <- mean(abs(fshd.best.orm.1pl.ld[lower.tri(fshd.best.orm.1pl.ld)]))
