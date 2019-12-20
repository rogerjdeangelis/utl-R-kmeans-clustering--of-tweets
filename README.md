# utl-R-kmeans-clustering--of-tweets
R kmeans clustering of tweets
    SAS and R kmeans clustering of tweets

    graphic output;
    https://tinyurl.com/r2ajxo4
    https://github.com/rogerjdeangelis/utl-R-kmeans-clustering--of-tweets/blob/master/kmeanscluster.pdf

    github
    https://github.com/rogerjdeangelis/utl-R-kmeans-clustering--of-tweets

    Python is probably a better language for this kind of problem, it scales better?

    https://cai.tools.sap/blog/text-clustering-with-r-an-introduction-for-data-scientists/

    *_                   _
    (_)_ __  _ __  _   _| |_
    | | '_ \| '_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    ;

    * words around keywords(ACA, FDA, Medicare and Cancer) in tweets;
    options validvarname=upcase;
    libname sd1 "d:/sd1";
    data sd1.have;
      infile cards delimiter=',';
      length tweets $44;
      informat tweets $44.;
      input tweets @@;
    cards4;
    health Medicare laws,first Medicare laws,health Medicare health,first Medicare laws
    health FDA laws,has FDA health,health FDA has,laws ACA health,first ACA has,the Cancer laws
    first Cancer has,health Cancer health,the Medicare the,the Medicare Law,Read Medicare Health
    Read FDA Law,Read FDA of,of ACA Read,Health ACA of,the ACA the,of ACA Law,the Cancer Read
    Read Cancer of,A Medicare of,of Medicare ACA,ACA Medicare of,of Medicare the,ACA FDA A
    preexisting FDA ACA,preexisting FDA preexisting,the ACA preexisting,ACA ACA ACA
    of ACA preexisting,the ACA A,A Cancer the,A Cancer A,ACA Cancer the,Young Medicare This
    of Medicare of,of Medicare ACA,Young FDA of,ACA FDA of,the ACA the,Young ACA the,of ACA the
    the Cancer the,This Cancer This,millions Medicare coverage,coverage Medicare to
    millions Medicare coverage,millions FDA expanded,millions FDA coverage,has FDA has
    coverage FDA expanded,has FDA has,millions ACA has,millions ACA has,coverage ACA expanded
    has ACA coverage,millions Cancer to,Lift Medicare on,on FDA on,Ban FDA Lift,on FDA Blood
    Ban FDA Ban,to ACA Ban,Lift Cancer Ban,to Cancer Lift,the Medicare On,the Medicare On
    Enrollment FDA On,Enrollment FDA the,Upswing FDA Upswing,Enrollment FDA Enrollment
    Enrollment FDA On,Enrollment Cancer Enrollment,the Cancer Upswing,FDA Medicare its
    decadeslong FDA ends,decadeslong FDA ban,ends FDA FDA,ban FDA decadeslong,decadeslong ACA its
    decadeslong ACA ban,decadeslong ACA FDA,ends ACA FDA,FDA Cancer its,ban Cancer its
    ban Cancer decadeslong,ban Cancer its,decadeslong Cancer ban,cancer Medicare cancer
    which Medicare skin,drug FDA which,which FDA cancer,works FDA skin,skin FDA which
    works ACA works,works ACA drug,skin Cancer skin,skin Cancer works,SCOTUS Medicare case
    case Medicare case,SCOTUS FDA v,Burwell FDA v,SCOTUS ACA on,Burwell ACA Burwell,v Cancer SCOTUS
    Burwell Cancer v,v Cancer on,CNBC Medicare too,The Medicare CNBC,too Medicare watches
    too Medicare CNBC,watches Medicare too,watches FDA too,watches FDA watches,The ACA CNBC
    CNBC ACA watches,The Cancer The,The Cancer The,CNBC Cancer backstory,too Cancer too
    watches Cancer watches,Internet Medicare Internet,the Medicare the,in Medicare Internet
    the Medicare Age,Age FDA in,Oversight FDA Internet,Internet FDA in,the FDA Oversight
    Oversight FDA Oversight,the ACA Internet,in Cancer Age,Age Cancer Age,Age Cancer the
    in Cancer Age,Approves Medicare Approves,Cancer Medicare AstraZeneca,Cancer Medicare Ovarian
    AstraZeneca Medicare AstraZeneca,AstraZeneca Medicare Ovarian,Ovarian FDA Drug
    Ovarian FDA Approves,Drug ACA Approves,Ovarian Cancer Drug,Cancer Cancer AstraZeneca
    Approves Cancer Approves,AstraZeneca Cancer AstraZeneca,Cancer Cancer Cancer
    Hepatitis Medicare Hepatitis,Backs Medicare Treatment,Hepatitis FDA Backs,AbbVies FDA Hepatitis
    Backs FDA Treatment,Backs Cancer Treatment,AbbVies Cancer Backs,Hepatitis Cancer Backs
    wins Medicare approval,wins Medicare FDA,wins FDA wins,approval FDA FDA,1st FDA FDA,wins FDA FDA
    wins FDA for,1st ACA wins,for Cancer FDA,1st Cancer approval,light Medicare green
    ABBV Medicare FDA,gets FDA FDA,light FDA green,light ACA green,ABBV ACA gets,light ACA green
    ABBV ACA light,green ACA FDA,light Cancer ABBV,green Cancer green,green Cancer green
    ABBV Cancer gets,light Cancer green,Humira Medicare Humira,awaits FDA Humira,built FDA Humira
    house FDA the,built FDA awaits,built FDA built,built ACA house,awaits ACA awaits,built ACA house
    awaits ACA house,built Cancer house,house Cancer house,awaits Cancer built,awaits Cancer the
    built Cancer built,on Medicare AbbVie,game Medicare AbbVie,pins Medicare hopes,game Medicare on
    pins FDA game,pins FDA AbbVie,AbbVie FDA pins,pins FDA pins,pins ACA on,AbbVie ACA game
    on ACA AbbVie,game Cancer game,game Cancer hopes,AbbVie Cancer game,game Cancer game
    on Medicare health,on FDA agency,agency FDA health,warns ACA health,warns Cancer warns
    warns Cancer on,agency Cancer top,top Cancer top,warns Cancer warns,Pharmalittle Medicare amp
    Good Medicare Good,Good Medicare Pharmalittle,Good Medicare Good
    Pharmalittle Medicare Pharmalittle,amp FDA amp,amp FDA Morning,Good FDA amp
    Good FDA Pharmalittle,Pharmalittle ACA sunshine,Good ACA amp,Morning ACA Pharmalittle
    Pharmalittle Cancer Pharmalittle,pure Medicare clear,says FDA pure,clear FDA clear,of FDA pure
    steer ACA steer,says ACA says,clear ACA clear,pure ACA says,pure ACA pure,of Cancer clear
    clear Cancer says,says Cancer steer,steer Cancer pure,steer Cancer steer,got Medicare FDA
    ok Medicare ok,FDA Medicare fast,scant FDA FDA,FDA FDA got,ok ACA got,scant ACA fast
    scant ACA got,fast ACA fast,got Cancer fast,got Cancer FDA,FDA Cancer scant,around Medicare 2
    morcellators Medicare morcellators,2 FDA morcellators,for FDA morcellators,for FDA 2,2 ACA 2
    morcellators Cancer were,around Cancer were,morcellators Cancer 2,for Cancer around
    were Cancer morcellators,Years Medicare did,Years Medicare did,19 Medicare did
    Years Medicare Years,19 FDA it,it FDA it,Years ACA Years,did Cancer 19,did Cancer did
    19 Cancer it,all Medicare doesnt,FDA Medicare FDA,all Medicare FDA,disclose Medicare FDA
    all Medicare financial,financial FDA FDA,disclose ACA financial,disclose ACA disclose
    FDA ACA doesnt,FDA ACA FDA,all Cancer all,FDA Cancer all,disclose Cancer FDA,all Cancer doesnt
    Mercks Medicare New,New Medicare Clears,Mercks Medicare Mercks,Mercks Medicare Mercks
    Clears FDA Mercks,of FDA of,Clears FDA Version,of ACA Mercks,Mercks ACA of,Version Cancer Clears
    that Medicare new,vaccine FDA Merck,vaccine FDA vaccine,vaccine FDA vaccine,new FDA vaccine
    clears FDA vaccine,new ACA Merck,vaccine ACA that,vaccine ACA that,that ACA that
    new Cancer Merck,vaccine Cancer Merck,panels Medicare panels,panels Medicare Some
    on FDA advising,advising FDA panels,on FDA advising,panels FDA panels,advising ACA docs
    on ACA Some,on ACA docs,docs ACA docs,on ACA on,on Cancer panels,panels Cancer docs
    13 Medicare of,panelists Medicare 13,on FDA of,13 FDA 13,panelists ACA FDA,FDA ACA FDA,13 ACA of
    FDA ACA of,13 ACA 13,13 Cancer on,panelists Cancer 13,docs Medicare FDA,advising FDA FDA
    advising FDA advising,of FDA of,advising FDA 10,advising ACA 10,10 ACA of,docs ACA FDA,10 ACA of
    10 Cancer 10,of Cancer advising,get Medicare you,on Medicare get,how Medicare on,do Medicare you
    you FDA do,get ACA you,get ACA do,on ACA you,on ACA you,you ACA on,you Cancer get,do Cancer how
    on Cancer on,on Cancer on,you Cancer do,FDA Medicare panels,judging Medicare panels
    FDA FDA devices,devices FDA devices,devices ACA FDA,judging ACA panels,for ACA FDA
    for Cancer FDA,advisory Medicare FDA,panels Medicare on,advisory FDA FDA,panels FDA panels
    doctors FDA advisory,doctors FDA FDA,on ACA FDA,doctors ACA panels,on ACA advisory
    panels ACA panels,advisory ACA FDA,doctors Cancer doctors,on Cancer on,FDA Cancer panels
    on Cancer on,on Cancer on,Financial Medicare Financial,Disclosed Medicare Disclosed
    Financial Medicare Ties,Not Medicare Advisers,Advisers Medicare Ties,Not FDA Disclosed
    Not FDA Financial,Ties FDA Disclosed,Not FDA Disclosed,Financial FDA Financial,Disclosed ACA Not
    Ties Cancer Not,Financial Cancer Ties,Not Cancer Ties,Financial Cancer Disclosed
    Disclosed Cancer Advisers
    ;;;;
    run;quit;



    SD1.HAVE total obs=408

    Obs    TWEETS

      1    health Medicare laws
      2    first Medicare laws
      3    health Medicare health
      4    first Medicare laws
      5    health FDA laws
      6    has FDA health
    ....


    *            _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| '_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    ;

             -0.25  -0.20  -0.15  -0.10  -0.05  0.00   0.05   0.10   0.15   0.20   0.25   0.30
            ---+------+------+------+------+------+------+------+------+------+------+------+---
         V1 |                                    . .                                           |
            |                                    . .                                    . .  . |
       0.20 +                                                                          . ..  . +  0.20
            |                    . .                                                           |
            |          . .                                                          . .. .     |
            |             .                                                          . .       |
            |       . .   .  . . .                                                             |
       0.15 +          . ...... .     . .                                                      +  0.15
            |             . . ...                                                . .           |
            |            . .  . . .                                                            |
            |                . ACA  . .                               . .                      |
            |             . .. . .... .. .                         . .  . ..   . .             |
       0.10 +                 . .. . . .                . .  . . .  . .  .. . . .              +  0.10
            |               . .. . ... .                    . . . . .....  . .                 |
            |                     . .  .                 . ..  . FDA . .. . . .                |
            |                      . . .. .              . . .  ......  . .                    |
            |                       . .. .               . .  ... .... . . . .                 |
       0.05 +                        . ... .                 . . ...... .                      +  0.05
            |                                               . .. .  . .                        |
            |                                                  . .                             |
            |                                               . .                                |
            |                      . .. .                                                      |
       0.00 +                    . .                                                           +  0.00
            |                           . . . .                                                |
            |                      . .  ......  . .     Note ACA and Medicare                  |
            |                   . .  ..... ...... .     are close together                     |
            |                         . Medicare .                                             |
      -0.05 +                         . ... . ... .                                            + -0.05
            |          . .             . .   . .                                               |
            |                               . .                                                |
            |                                                                                  |
            |                                                             . .. . .             |
      -0.10 +                                                                 . .   . .        + -0.10
            |                                                                 . .              |
            |                                       . .. . .                   . .   . .       |
            |                                     . .  . .                                     |
            |                                       . ... .                                    |
      -0.15 +                                       . . ... .                                  + -0.15
            |                                      . . . .. .                                  |
            |                                     . . .....  . .                               |
            |                                        .CANCER. .                                |
            |                                         . .  ... .                               |
      -0.20 +                                         . ... .                                  + -0.20
            |                                            . .     . .                           |
            |                                       ... .                                      |
            |                                  . .   . . .. .                                  |
            |             . .                          . .  . .                                |
      -0.25 +                                 . .                                              + -0.25
            |                                           . .                                    |
            ---+------+------+------+------+------+------+------+------+------+------+------+---
             -0.25  -0.20  -0.15  -0.10  -0.05  0.00   0.05   0.10   0.15   0.20   0.25   0.30



    *          _       _   _
     ___  ___ | |_   _| |_(_) ___  _ __
    / __|/ _ \| | | | | __| |/ _ \| '_ \
    \__ \ (_) | | |_| | |_| | (_) | | | |
    |___/\___/|_|\__,_|\__|_|\___/|_| |_|

    ;

    * ensure a rerun does not use previous tables;

    proc datasets lib=work nolist;
    delete want wntfin wntCmb;
    run;quit;

    %utlfkil(d:/xpt/want.xpt);

    %utl_submit_r64('
    library(tm);
    library(haven);
    library(proxy);
    library(dbscan);
    library(data.table);
    library(SASxport);
    have<-read_sas("d:/sd1/have.sas7bdat");
    truth.K<-4;
    corpus.cleaned = tm::Corpus(tm::VectorSource(have$TWEETS));
    tdm <- tm::DocumentTermMatrix(corpus.cleaned);
    tdm.tfidf <- tm::weightTfIdf(tdm);
    tdm.tfidf <- tm::removeSparseTerms(tdm.tfidf, 0.999);
    tfidf.matrix <- as.matrix(tdm.tfidf);
    dist.matrix = proxy::dist(tfidf.matrix, method = "cosine");
    clustering.kmeans <- kmeans(tfidf.matrix, 3);
    clustering.hierarchical <- hclust(dist.matrix, method = "ward.D2");
    clustering.dbscan <- dbscan::hdbscan(dist.matrix, minPts = 10);
    master.cluster <- clustering.kmeans$cluster;
    slave.hierarchical <- cutree(clustering.hierarchical, k = truth.K);
    slave.dbscan <- clustering.dbscan$cluster;
    stacked.clustering <- rep(NA, length(master.cluster));
    names(stacked.clustering) <- 1:length(master.cluster);
    for (cluster in unique(master.cluster)) {
      indexes = which(master.cluster == cluster, arr.ind = TRUE);
      slave1.votes <- table(slave.hierarchical[indexes]);
      slave1.maxcount <- names(slave1.votes)[which.max(slave1.votes)];
      slave1.indexes = which(slave.hierarchical == slave1.maxcount, arr.ind = TRUE);
      slave2.votes <- table(slave.dbscan[indexes]);
      slave2.maxcount <- names(slave2.votes)[which.max(slave2.votes)];
      stacked.clustering[indexes] <- slave2.maxcount;
    };
    points <- cmdscale(dist.matrix, k = truth.K);
    want<-as.data.table(points);
    previous.par <- par(mfrow=c(2,2), mar = rep(1.5, 4));
    pdf(file="d:/pdf/kmeanscluster.pdf");
    palette <- colorspace::diverge_hcl(truth.K);
    plot(points, main = "K-Means clustering", col = as.factor(master.cluster),
         mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0));
    plot(points, main = "Hierarchical clustering", col = as.factor(slave.hierarchical),
         mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0));
    plot(points, main = "Density-based clustering", col = as.factor(slave.dbscan),
         mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0));
    plot(points, main = "Stacked clustering", col = as.factor(stacked.clustering),
         mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0));
    par(previous.par);
    pdf();
    write.xport(want,file="d:/xpt/want.xpt");
    ');

    libname xpt xport "d:/xpt/want.xpt";
    data want;
      retain label " ";
      merge sd1.have xpt.want;
      keyword=scan(tweets,2);
      keep KEYWORD V1 V2 V3 V4;
    run;quit;
    libname xpt clear;

    proc summary data=want mean nway;
    class keyword;
    var v1 v2 v3 v4;
    output out=wntfin(DROP=_:) mean=;
    run;quit;

    data wntCmb;
       length label $3;
       set want wntfin(in=b);
       if b then label=keyword;
       else label='*';
    run;quit;

    options ls=96 ps=64;
    proc plot data=wntCmb;
       title "K-Means clustering";
       plot v1*v2='*' $ label /haxis =-.25 to .3 by .05 vaxis =-.25 to .2 by .05  box ;
    run;quit;


