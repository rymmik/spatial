library(sf)
library(ggplot2)
library(sp)
library(spatialreg)
library(spdep)
library(RColorBrewer)
library(classInt)
library(cartography)
library(tmap)
library(osmdata)
library(spgwr)
library(factoextra)
library(readxl)
library(GWmodel)
library(terra)	
library(fossil)	
library(plot.matrix)
library(doBy)

options(scipen = 999)

#Dane do modelu
dane<-read.csv("D:/Studia/WNE/Materialy/Rok 4/Ekonometria przestrzenna/projekt/dane_nowe.csv")

#Dane shp
setwd("D:/Studia/WNE/Materialy/Rok 4/Ekonometria przestrzenna/dane")
POW<-st_read("dzielnice_Warszawy/dzielnice_Warszawy.shp")
POW<-st_transform(POW, 4326) 	# konwersja 4326=WGS84, 4267=NAD27 
pow<-as_Spatial(POW, cast = TRUE, IDs ="jpt_kod_je")
pow<-spTransform(pow, CRS("+proj=longlat +datum=NAD83"))

#Mapa cen
cols<-rainbow(4)
stat<-summary(dane$cena_m2)
brks<-c(stat[1], stat[2],stat[4],stat[5], stat[6])
plot(st_geometry(POW))
points(dane[,c("wspolrzedne_wysokosc", "wspolrzedne_szerokosc")], pch=21, bg=cols[findInterval(dane$cena_m2, brks)])
legend("topright", legend = c(paste("<", round(stat[2], 0)),
                              paste(round(stat[2],0), "to", round(stat[4], 0)),
                              paste(round(stat[4],0), "to", round(stat[5], 0)),
                              paste(">", round(stat[5], 3))),
       fill = cols[1:4], title = paste("Cena za m2"))

#przeksztalcenia
dane$wiek_budynku<-2023-dane$rok_budowy
dane$cena_m2<-log(dane$cena_m2)

# losowe sortowanie zbioru danych
dane$los<-runif(dim(dane)[1], 0,1)
dane<-orderBy(~los, data=dane)

# losowa korekta zmiennych (jitter) (bo wartości są zbyt podobne)
dane$los2<-rnorm(dim(dane)[1], 0,0.0005)
dane$wspolrzedne_szerokosc<-dane$wspolrzedne_szerokosc+dane$los2
dane$los3<-rnorm(dim(dane)[1], 0,0.0005)
dane$wspolrzedne_wysokosc<-dane$wspolrzedne_wysokosc+dane$los3


#Macierz sąsiedztwa dla sąsiadów w promieniu d km
crds<-cbind(dane$wspolrzedne_wysokosc, dane$wspolrzedne_szerokosc)	# współrzędne centroidów / centroids
conti.nb<-dnearneigh(crds, 0, 3,longlat=TRUE)	# w km / in km
conti.listw<-nb2listw(conti.nb, zero.policy=TRUE)
conti.m<-nb2mat(conti.nb, zero.policy=TRUE)
conti.m[1:10, 1:10]
summary(conti.nb)
#K sasaidow
conti<-knearneigh(crds, k=5)
conti.nb<-knn2nb(conti)
conti.listw<-nb2listw(make.sym.nb(conti.nb))
conti.m<-nb2mat(conti.nb, zero.policy=TRUE)
conti.m[1:10, 1:10]
summary(conti.nb)



# I Morana dla macierzy odwrotnej odległości
# Moran’s I for inverse distance matrix
result<-moran.test(dane$cena_m2, conti.listw)
result

# wykres punktowy Morana – krok po kroku
# Moran scatterplot – step by step version
x<-dane$cena_m2 # wybór zmiennej / variable selection
zx<-as.data.frame(scale(x))  # standaryzacja / standardization of variable
wzx<-lag.listw(conti.listw, zx$V1) #spatial lag of x / opóźnienie przestrzenne x
morlm<-lm(wzx~zx$V1) # linear regression  / regresja liniowa           
summary(morlm)
slope<-morlm$coefficients[2] # coefficient from regression / współczynnik kier.
intercept<-morlm$coefficients[1] # constant term in regression / stała z regr.
plot(zx$V1, wzx, xlab="zx",ylab="spatial lag of zx", pch="*")  
abline(intercept, slope)  # regression line / linia regresji
abline(h=0, lty=2) # supplementary horizontal line y=0 / linia pozioma
abline(v=0, lty=2) # supplementary vertical line x=0 / linia pionowa

#Join-count test 
var.factor<-factor(cut(dane$cena_m2, breaks=c(0,10000, 17000,100000), labels=c("low", "medium", "high"))) # factor–values into labels / factor-wartości na etykiety
head(var.factor)
brks1<-c(0,10000, 17000,100000)
cols<-c("green", "blue", "red")
plot(dane$cena_m2, bg=cols[findInterval(dane$cena_m2, brks1)], pch=21)
abline(h=c(0,10000, 17000,100000), lty=3)
joincount.test(var.factor, conti.listw)


#Mapa punktów
plot(st_geometry(POW))
points(dane[, c("wspolrzedne_wysokosc", "wspolrzedne_szerokosc")], pch=".", cex=1.5)





# modele
eq<-cena_m2~remont + wykonczenie + winda + balkon + taras + wiek_budynku + spoldzielcze + powierzchnia  + pietro + pietro_max  + odleglosc_metro + odleglosc_centrum + agencja

#MNK
model.lm<-lm(eq, data=dane)
summary(model.lm)

#GWR
#optymalne pasmo
daneGWR<-dane
dMat=gw.dist(crds)
coordinates(daneGWR)<-c("wspolrzedne_wysokosc", "wspolrzedne_szerokosc")
bw<-GWmodel::bw.gwr(eq, data=daneGWR, kernel='gaussian', adaptive=FALSE, dMat=dMat)
#bw=0.3949375
#model
modelGWR<-GWmodel::gwr.basic(eq, data=daneGWR, kernel='gaussian', adaptive=FALSE, bw=0.3949375, dMat=dMat)  
modelGWR
head(modelGWR$SDF)
modelGWR$SDF[,2:14]

# wykresy pudełkowe (boxplots) współczynników – aby zobaczyc zróżnicowanie GWR
boxplot(as.data.frame(modelGWR$SDF)[,2:14])
abline(h=0, lty=3, lwd=2, col="red")


#Wykresy współczynnikow GWR dla zmiennej
names(modelGWR$SDF)
zmienna<-as.data.frame(modelGWR$SDF[,1])
cols<-rainbow(10)
stat<-summary(zmienna[,1])
brks<-c(stat[1], stat[4]-sd(zmienna[,1]),stat[4]+sd(zmienna[,1]),stat[6])
plot(st_geometry(POW))
points(dane[,c("wspolrzedne_wysokosc", "wspolrzedne_szerokosc")], pch=21, bg=cols[findInterval(zmienna[,1], brks)])

legend("topright", legend = c(paste("<", round(stat[4] - sd(zmienna[, 1]), 3)),
                              paste(round(stat[4] - sd(zmienna[, 1]), 3), "to", round(stat[4] + sd(zmienna[, 1]), 3)),
                              paste(">", round(stat[4] + sd(zmienna[, 1]), 3))),
       fill = cols[1:3], title = paste("GWR",names(zmienna[1])))

#Klastry
fviz_nbclust(as.data.frame(modelGWR$SDF[,2:14]), FUNcluster=kmeans, method = 'silhouette')
fviz_nbclust(as.data.frame(modelGWR$SDF[,2:14]), FUNcluster = kmeans, method = 'wss')
fviz_nbclust(as.data.frame(modelGWR$SDF[,2:14]), FUNcluster = kmeans, method = 'gap_stat')
klastry_b<-eclust(as.data.frame(modelGWR$SDF[,2:14]), "kmeans", k=8) #liczba klastrów
fviz_silhouette(klastry_b)
klastry_b<-kmeans(as.data.frame(modelGWR$SDF[,2:14]), 8) 
zmienna<-klastry_b$cluster
cols<-rainbow(8)
brks<-c(1,2,3,4,5,6,7,8)
plot(st_geometry(POW))
points(dane[,c("wspolrzedne_wysokosc", "wspolrzedne_szerokosc")], pch=21, bg=cols[findInterval(zmienna, brks)])



#Rastrowanie
bb<-st_bbox(POW)
r<-rast(nrows=50, ncols=50, xmin=bb[1], ymin=bb[2], xmax=bb[3], ymax=bb[4])
raster<-rasterize(as.matrix(crds), r, value=klastry_b$clust,  fun=median) 
plot(raster)
plot(st_geometry(POW), add=TRUE)
#Przypisanie klastrów do obserwacji
dane$klaster1<-rep(0, times=dim(dane)[1])
dane$klaster1[klastry_b$cluster==1]<-1
dane$klaster2<-rep(0, times=dim(dane)[1])
dane$klaster2[klastry_b$cluster==2]<-1
dane$klaster3<-rep(0, times=dim(dane)[1])
dane$klaster3[klastry_b$cluster==3]<-1
dane$klaster4<-rep(0, times=dim(dane)[1])
dane$klaster4[klastry_b$cluster==4]<-1
dane$klaster5<-rep(0, times=dim(dane)[1])
dane$klaster5[klastry_b$cluster==5]<-1
dane$klaster6<-rep(0, times=dim(dane)[1])
dane$klaster6[klastry_b$cluster==6]<-1
dane$klaster7<-rep(0, times=dim(dane)[1])
dane$klaster7[klastry_b$cluster==7]<-1


#Nowe modele z klastrami
eq<-cena_m2~remont + wykonczenie + winda + balkon + taras + wiek_budynku + spoldzielcze + powierzchnia  + pietro + pietro_max  + odleglosc_metro + odleglosc_centrum + agencja +  klaster1 + klaster2 + klaster3 + klaster4 + klaster5 + klaster6 + klaster7 

#Model liniowy
model.ols<-lm(eq, data=dane)
summary(model.ols)

#Model SEM
model.sem<-errorsarlm(eq, data=dane, conti.listw)
summary(model.sem)








