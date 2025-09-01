#include "ExprPredictor.h"
#include <fstream>
#include <sstream>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
string trim(const string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}
bool parseConfigFile(const string& configFile, map<string, string>& config) {
    ifstream file(configFile.c_str());
    if (!file.is_open()) {
        return false;
    }
    string line;
    while (getline(file, line)) {
        // Remove comments
        size_t commentPos = line.find('#');
        if (commentPos != string::npos) {
            line = line.substr(0, commentPos);
        }
        
        line = trim(line);
        if (line.empty()) continue;
        
        size_t pos = line.find('=');
        if (pos != string::npos) {
            string key = trim(line.substr(0, pos));
            string value = trim(line.substr(pos + 1));
            config[key] = value;
        }
    }
    file.close();
    return true;
}
int main( int argc, char* argv[] ) 
{
    string seqFile, annFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, seqFileb, bIFile, dcFile, duFile, adamiFile, exprFile2;
    string outFile;     
    double coopDistThr = 150;
    double factorIntSigma = 50.0;   
    double repressionDistThr = 150;
    double energyThr = 5;
    int maxContact = 1;
    string configFile;
    bool useConfig = false;
    for (int i = 1; i < argc; i++) {
        if (!strcmp("-config", argv[i]) || !strcmp("-conf", argv[i])) {
            if (i + 1 < argc) {
                configFile = argv[++i];
                useConfig = true;
                break;
            }
        }
    }
    map<string, string> config;
    if (useConfig) {
        if (!parseConfigFile(configFile, config)) {
            cerr << "Error: Cannot open config file " << configFile << endl;
            exit(1);
        }
        if (config.find("seqFile") != config.end()) seqFile = config["seqFile"];
        if (config.find("annFile") != config.end()) annFile = config["annFile"];
        if (config.find("exprFile") != config.end()) exprFile = config["exprFile"];
        if (config.find("motifFile") != config.end()) motifFile = config["motifFile"];
        if (config.find("factorExprFile") != config.end()) factorExprFile = config["factorExprFile"];
        if (config.find("coopFile") != config.end()) coopFile = config["coopFile"];
        if (config.find("factorInfoFile") != config.end()) factorInfoFile = config["factorInfoFile"];
        if (config.find("repressionFile") != config.end()) repressionFile = config["repressionFile"];
        if (config.find("paramFile") != config.end()) parFile = config["paramFile"];
        if (config.find("bindingSitesFile") != config.end()) seqFileb = config["bindingSitesFile"];
        if (config.find("bindingIntensityFile") != config.end()) bIFile = config["bindingIntensityFile"];
        if (config.find("dcFile") != config.end()) dcFile = config["dcFile"];
        if (config.find("duFile") != config.end()) duFile = config["duFile"];
        if (config.find("adamiFile") != config.end()) adamiFile = config["adamiFile"];
        if (config.find("exprFile2") != config.end()) exprFile2 = config["exprFile2"];
        if (config.find("outputFile") != config.end()) outFile = config["outputFile"];
        if (config.find("coopDistThr") != config.end()) coopDistThr = atof(config["coopDistThr"].c_str());
        if (config.find("factorIntSigma") != config.end()) factorIntSigma = atof(config["factorIntSigma"].c_str());
        if (config.find("repressionDistThr") != config.end()) repressionDistThr = atof(config["repressionDistThr"].c_str());
        if (config.find("energyThreshold") != config.end()) energyThr = atof(config["energyThreshold"].c_str());
        if (config.find("maxContact") != config.end()) maxContact = atoi(config["maxContact"].c_str());
        if (config.find("modelOption") != config.end()) ExprPredictor::modelOption = getModelOption(config["modelOption"].c_str());
        if (config.find("objOption") != config.end()) ExprPredictor::objOption = getObjOption(config["objOption"].c_str());
    }
    bool free_fix_indicators[] = {0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    0,0,0 };
    int nExps = 1;	
    vector <bool> indicator_bool ( free_fix_indicators, free_fix_indicators + sizeof( free_fix_indicators )/sizeof( bool ));
    int binwidth; 
    ExprPredictor::nAlternations = 1;
    int nalt =1;
    int nrand2=1;
    if (useConfig) {
        if (config.find("nAlternations") != config.end()) nalt = atoi(config["nAlternations"].c_str());
        if (config.find("nRandomStarts") != config.end()) nrand2 = atoi(config["nRandomStarts"].c_str());
        if (config.find("nExps") != config.end()) nExps = atoi(config["nExps"].c_str());
        if (config.find("binwidth") != config.end()) binwidth = atoi(config["binwidth"].c_str());
    }
    for ( int i = 1; i < argc; i++ ) {
        if (!strcmp("-config", argv[i]) || !strcmp("-conf", argv[i])) {
            i++; 
            continue;
        }
        if ( !strcmp( "-s", argv[ i ] ) )
        seqFile = argv[ ++i ];
        else if ( !strcmp( "-nrand", argv[ i ] ) )
        nrand2 = atoi( argv[++i] );
        else if ( !strcmp( "-a", argv[ i ] ) )
        annFile = argv[ ++i ];
        else if ( !strcmp( "-du", argv[ i ] ) )
        duFile = argv[ ++i ];  
        else if ( !strcmp( "-dc", argv[ i ] ) )
        dcFile = argv[ ++i ];
        else if ( !strcmp( "-sa", argv[ i ] ) )
        adamiFile = argv[ ++i ];            
        else if ( !strcmp( "-e", argv[ i ] ) )
        exprFile = argv[ ++i ];   
        else if ( !strcmp( "-e2", argv[ i ] ) )
        exprFile2 = argv[ ++i ];          
        else if ( !strcmp( "-m", argv[ i ] ) )
        motifFile = argv[ ++i ];
        else if ( !strcmp( "-f", argv[ i ] ) )
        factorExprFile = argv[ ++i ];    
        else if ( !strcmp( "-o", argv[ i ] ) )
        ExprPredictor::modelOption = getModelOption( argv[++i] );
        else if ( !strcmp( "-c", argv[ i ] ) )
        coopFile = argv[ ++i ];
        else if ( !strcmp( "-i", argv[ i ] ) )
        factorInfoFile = argv[ ++i ];            
        else if ( !strcmp( "-r", argv[ i ] ) )
        repressionFile = argv[ ++i ];  
        else if ( !strcmp( "-oo", argv[ i ] ) )
        ExprPredictor::objOption = getObjOption( argv[++i] );    
        else if ( !strcmp( "-mc", argv[i] ) )
        maxContact = atoi( argv[++i] );
        else if ( !strcmp( "-fo", argv[i] ) )
        outFile = argv[++i];    
        else if ( !strcmp( "-p", argv[i] ) )
        parFile = argv[++i]; 
        else if ( !strcmp( "-ct", argv[i] ) )
        coopDistThr = atof( argv[++i] ); 
        else if ( !strcmp( "-binwt", argv[i] ) )
        binwidth = atof( argv[++i] ); 
        else if ( !strcmp( "-sigma", argv[i] ) )
        factorIntSigma = atof( argv[++i] );    
        else if ( !strcmp( "-rt", argv[i] ) )
        repressionDistThr = atof( argv[++i] );    
        else if ( !strcmp( "-et", argv[i] ) )
        energyThr = atof( argv[++i] );
        else if ( !strcmp( "-na", argv[i] ) )
        nalt = atoi( argv[++i] );  
        else if ( !strcmp( "-n", argv[ i ] ) )
        nExps = atoi( argv[ ++i ] );
        else if ( !strcmp( "-bs", argv[i] ) )
        seqFileb = argv[++i];    
        else if ( !strcmp( "-bi", argv[i] ) )
        bIFile = argv[++i]; 
    }
    if ( seqFile.empty() || exprFile.empty() || motifFile.empty() || factorExprFile.empty() || outFile.empty() || ( ( ExprPredictor::modelOption == QUENCHING || ExprPredictor::modelOption == CHRMOD_UNLIMITED || ExprPredictor::modelOption == CHRMOD_LIMITED ) &&  factorInfoFile.empty() ) || ( ExprPredictor::modelOption == QUENCHING && repressionFile.empty() ) ) {
        cerr << "Usage: " << argv[ 0 ] << " -config configFile" << endl;
        cerr << "   OR: " << argv[ 0 ] << " -s seqFile -e exprFile -m motifFile -f factorExprFile -fo outFile [-a annFile -o modelOption -c coopFile -i factorInfoFile -r repressionFile -oo objOption -mc maxContact -p parFile -rt repressionDistThr -et energyThr -na nAlternations -ct coopDistThr -sigma factorIntSigma]" << endl;
        cerr << "modelOption: Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited" << endl;
        cerr << "objOption: corr, SSE, Cross_Corr" << endl;
        exit( 1 );
    }           
    vector< int > mmm(4);
    mmm[0] = 30;
    mmm[1] = 52;
    mmm[2] = 155;
    mmm[3] = 160;
    vector< int > mmmr(4);
    mmmr[0] = 50;
    mmmr[1] = 100;
    mmmr[2] = 150;
    mmmr[3] = 200;
    double gcContent = 0.4;
    FactorIntType intOption = BINSF;     
    ExprPar::searchOption =  CONSTRAINED;      
    ExprPar::estBindingOption = 1;
    int nbin = mmm.size() + 1;
    ExprPar::nbins = nbin;    
    ExprPredictor::nAlternations=1;
    ExprPredictor::nRandStarts = 1;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_CrossCorr = 1.0E-10;
    ExprPredictor::nSimplexIters = 3;
    ExprPredictor::nGradientIters = 10;
    
    // Create oData directory if it doesn't exist
    struct stat st = {0};
    if (stat("oData", &st) == -1) {
        mkdir("oData", 0700);
    }
    
    int rval;
    vector< vector< double > > data;    
    vector< string > labels;    
    string factor1, factor2;    
    vector< Sequence > seqs;
    vector< string > seqNames;
    rval = readSequences( seqFile, seqs, seqNames );
    assert( rval != RET_ERROR );
    int nSeqs = seqs.size();
    rval = readSequences( seqFile, ExprPredictor::seqsy, ExprPredictor::seqNmes );
    vector< string > condNames;  
    rval = readMatrix( exprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nSeqs );
    for ( int i = 0; i < nSeqs; i++ ) assert( labels[i] == seqNames[i] );
    Matrix exprData( data ); 
    int nConds = exprData.nCols();
    ExprPredictor::exprData2 = exprData;
    vector< Motif > motifs;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );
    rval = readMotifs( motifFile, background, motifs, motifNames ); 
    assert( rval != RET_ERROR );
    int nFactors = motifs.size();
    Motif d;
    d.copy(motifs[0]);
    Matrix dtemp=d.getPwm2();
    vector< double > e(4);
    e[0] = 1;
    e[1] = 1;
    e[2] = 1;
    e[3] = 10;
    vector< vector< double > > d2s;
    for(int i = 0 ; i < dtemp.nRows(); i++ ) {
        d2s.push_back( dtemp.getRow(i) ) ;
        if ( i == 4 ){
            d2s.push_back( e );
        }
    }
    Matrix d2sm( d2s );
    Matrix dcm =countmatrixS( dcFile );
    Motif dcmm( dcm, background);
    Matrix dum = countmatrixS( duFile );
    Motif d2( dum, background);
    map< string, int > factorIdxMap;
    for ( int i = 0; i < motifNames.size(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
    labels.clear();
    data.clear();
    rval = readMatrix( factorExprFile, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nFactors && condNames.size() == nConds );
    for ( int i = 0; i < nFactors; i++ ) assert( labels[i] == motifNames[i] );
    Matrix factorExprData( data ); 
    assert( factorExprData.nCols() == nConds ); 
    ExprPredictor::factorExprData2 = factorExprData;
    vector< double > energyThrs( nFactors, energyThr ); 
    vector< SiteVec > seqSites( nSeqs );
    vector< int > seqLengths( nSeqs );
    SeqAnnotator ann( motifs, energyThrs );
    if ( annFile.empty() ) {        
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.annot( seqs[ i ], seqSites[ i ] );
            seqLengths[i] = seqs[i].size();
        }
    } else {    
        rval = readSites( annFile, factorIdxMap, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }
    vector< SiteVec > seqSitesbot( nSeqs );
    seqSitesbot = seqSites;
    vector< SiteVec > seqSitesm1( nSeqs );
    seqSitesm1 = seqSites;
    vector< SiteVec > seqSitesm2( nSeqs );
    seqSitesm2 = seqSites;
    vector< SiteVec > seqSitesf2( nSeqs );
    seqSitesf2 = seqSites;
    vector< SiteVec > seqSitesbotf2(nSeqs);
    seqSitesbotf2 = seqSites;
    vector< SiteVec > seqSitesm1f2(nSeqs);
    seqSitesm1f2 = seqSites; 
    vector< SiteVec > seqSitesm2f2( nSeqs);
    seqSitesm2f2 = seqSites;
    vector< SiteVec >  seqSitesf3( nSeqs);
    seqSitesf3 = seqSites;
    vector< SiteVec >  seqSitesbotf3( nSeqs);
    seqSitesbotf3 = seqSites;
    vector< SiteVec > seqSitesm1f3( nSeqs );
    seqSitesm1f3 = seqSites;
    vector< SiteVec >  seqSitesm2f3( nSeqs);
    seqSitesm2f3=seqSites;
    vector< vector< Site > > seqSitesm1delete1;
    seqSitesm1delete1 = seqSites;
    SeqAnnotator anny( motifs, energyThrs , d2);
    IntMatrix coopMat( nFactors, nFactors, false );
    if ( !coopFile.empty() ) {
        ifstream fcoop( coopFile.c_str() );
        if ( !fcoop ) {
            cerr << "Cannot open the cooperativity file " << coopFile << endl;
            exit( 1 );
        }  
        while ( fcoop >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            coopMat( idx1, idx2 ) = true;
            coopMat( idx2, idx1 ) = true;
        }        
    } 
    vector< bool > actIndicators( nFactors, false );
    vector< bool > repIndicators( nFactors, false );
    if ( !factorInfoFile.empty() ) {
        ifstream finfo( factorInfoFile.c_str() );
        if ( !finfo ) {
            cerr << "Cannot open the factor information file " << factorInfoFile << endl;
            exit( 1 );
        }      
        string name;
        int i = 0, actRole, repRole;
        while ( finfo >> name >> actRole >> repRole ) {
            assert( name == motifNames[i] );
            if ( actRole ) actIndicators[i] = true;
            if ( repRole ) repIndicators[i] = true;
            i++;
        }
    }
    IntMatrix repressionMat( nFactors, nFactors, false );
    if ( !repressionFile.empty() ) {
        ifstream frepr( repressionFile.c_str() );
        if ( !frepr ) {
            cerr << "Cannot open the repression file " << repressionFile << endl;
            exit( 1 );
        }        
        while ( frepr >> factor1 >> factor2 ) {
            assert( factorIdxMap.count( factor1 ) && factorIdxMap.count( factor2 ) );
            int idx1 = factorIdxMap[factor1];
            int idx2 = factorIdxMap[factor2];
            repressionMat( idx1, idx2 ) = true;
        }        
    }
    vector< Sequence > allSeqs;
    vector< string > allNames;
    rval = readSequences( seqFileb, allSeqs, allNames );
    assert( rval != RET_ERROR );
    int nSeqsb = allSeqs.size() / nExps;
    vector< vector< Sequence > > seqs2( nExps );
    vector< vector< string > > names( nExps );	
    int counter = 0;
    for ( int i = 0; i < nExps; i++ ) {
        for ( int j = 0; j < nSeqsb; j++ ) {
            seqs2[ i ].push_back( allSeqs[ counter ] );	
            names[ i ].push_back( allNames[ counter ] );
            counter++;
        }
    }
    ifstream fdata( bIFile.c_str() );
    if ( !fdata ) {
        cerr << "Cannot find the binding data file " << bIFile << endl;
        exit( 1 );
    }
    vector< vector< double > > bindingData( nExps, vector<double>(nSeqsb) );
    for ( int i = 0; i < nExps; i++ ) {
        for ( int j = 0; j < nSeqsb; j++ ) {
            string name;
            fdata >> name;			
            if( name != names[ i ][ j ] ) { 
                cerr << "Error: " << names[ i ][ j ] << "\t" << name << endl;
                exit( 1 );
            }
            fdata >> bindingData[ i ][ j ];	
        }
    }	
    vector< vector< SiteVec > > seqSitesb( seqs2.size() );  
    SeqAnnotator annb( motifs, energyThrs );
    for ( int i = 0; i < nExps; i++ ) {
        for ( int j = 0; j < nSeqsb; j++ ) {
            SiteVec sites;
            annb.annot( seqs2[ i ][ j ], sites );  
            seqSitesb[ i ].push_back( sites );
        }
    }	
    FactorIntFunc* intFunc; 
    if ( intOption == BINARY ) intFunc = new FactorIntFuncBinary( coopDistThr ); 
    else if ( intOption == GAUSSIAN ) intFunc = new FactorIntFuncGaussian( coopDistThr, factorIntSigma );
    else if ( intOption == BINSF ) intFunc = new FactorIntFuncBinsf( coopDistThr, nbin );
    else {
        cerr << "Interaction Function is invalid " << endl; exit( 1 ); 
    }
    ExprPredictor* predictor = new ExprPredictor(seqSitesb, bindingData, seqSites, seqLengths,
    exprData, motifs, factorExprData, intFunc, coopMat, actIndicators, maxContact, repIndicators,
    repressionMat, repressionDistThr, mmm, mmmr, indicator_bool, anny, exprFile, seqSitesbot, 
    seqSitesm1,seqSitesm2, seqSitesf2 ,seqSitesbotf2, seqSitesm1f2 ,seqSitesm2f2, seqSitesf3, 
    seqSitesbotf3,seqSitesm1f3, seqSitesm2f3 );  
    ExprPar par_init( nFactors);  
    if ( !parFile.empty() ) {
        rval = par_init.load( parFile, coopMat, repIndicators );
        if ( rval == RET_ERROR ) {
            cerr << "Cannot read parameters from " << parFile << endl;
            exit( 1 );
        } 
    }
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;	
    rng = gsl_rng_alloc( T );
    gsl_rng_set( rng, 12345 );		
    (*predictor).clasvar = 1;
cout<< "out" << endl;
    predictor->objFuncborder(  par_init );
    double obj = predictor->getObj();
    predictor->printPar(predictor->getPar());
    (*predictor).clasvar = 0;
    ExprPar::searchOption =  UNCONSTRAINED;
    bool free_fix_indicators2[] = {1,0,1,
    0,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    1,0,0 };
    vector <bool> indicator_bool2 ( free_fix_indicators2, free_fix_indicators2 + sizeof( free_fix_indicators2 )/sizeof( bool ));
    predictor->modifyIndicatorbool( indicator_bool2);
    ExprPar::searchOption =  UNCONSTRAINED;
    ExprPredictor::nRandStarts = nrand2;
    ExprPredictor::nAlternations=nalt;
    ExprPredictor::nSimplexIters = 20;
    ExprPredictor::nGradientIters = 10;
    predictor->printPar(predictor->getPar());
    double obj_before = predictor->getObj();
    predictor->train( predictor->getPar(), rng);
    // predictor->train4( predictor->getPar() );
    predictor->printPar(predictor->getPar());
    obj = predictor->getObj();
    double improvement = abs(obj - obj_before);
    if (improvement < 1e-6) {
    }
    ExprPar::searchOption =  UNCONSTRAINED;
    ofstream to3("oData/ot3.txt");
    to3.close();
    predictor->printPar(predictor->getPar());
    ofstream fo("oData/pars2.txt");
    fo.close();
    seqs.clear();
    seqNames.clear();
    bool free_fix_indicators3[] = {0,0,0,
    1,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    1,0,1 };
    vector <bool> indicator_bool3 ( free_fix_indicators3, free_fix_indicators3 + sizeof( free_fix_indicators3 )/sizeof( bool ));
    if (adamiFile.empty()) {
        adamiFile = seqFile;
    }
    rval = readSequences(adamiFile, seqs, seqNames );
    assert( rval != RET_ERROR );
    nSeqs = seqs.size();
    data.clear();
    labels.clear(); 
    if (exprFile2.empty()) {
        exprFile2 = exprFile;
    }
    rval = readMatrix( exprFile2, labels, condNames, data );
    assert( rval != RET_ERROR );
    assert( labels.size() == nSeqs );
    for ( int i = 0; i < nSeqs; i++ ) assert( labels[i] == seqNames[i] );
    vector< vector< double > > data2;    
    vector< string > labels2; 
    Matrix exprData2( data ); 
    nConds = exprData2.nCols();
    seqSites.clear();
    vector< SiteVec > seqSitestes( nSeqs );
    seqSites=seqSitestes;
    vector< int > seqLengthstes( nSeqs );
    seqLengths.clear();
    seqLengths=seqLengthstes;
    if ( annFile.empty() ) {        
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.annot( seqs[ i ], seqSites[ i ] );
            seqLengths[i] = seqs[i].size();
        }
    } else {    
        rval = readSites( annFile, factorIdxMap, seqSites, true );
        assert( rval != RET_ERROR );
        for ( int i = 0; i < nSeqs; i++ ) {
            ann.compEnergy( seqs[i], seqSites[i] );
            seqLengths[i] = seqs[i].size();
        }
    }



    seqSitesbot.clear(); 
    seqSitesm1.clear(); 
    seqSitesm2.clear();
    seqSitesf2.clear();
    seqSitesbotf2.clear();
    seqSitesm1f2.clear();
    seqSitesm2f2.clear();
    seqSitesf3.clear();
    seqSitesbotf3.clear();
    seqSitesm1f3.clear();
    seqSitesm2f3.clear();
    seqSitesm1delete1.clear();
    seqSitesbot=seqSites;
    seqSitesm1 = seqSites;
    seqSitesm2 = seqSites;
    seqSitesf2 = seqSites;
    seqSitesbotf2 = seqSites;
    seqSitesm1f2 = seqSites; 
    seqSitesm2f2 = seqSites;
    seqSitesf3 = seqSites;
    seqSitesbotf3 = seqSites;
    seqSitesm1f3 = seqSites;
    seqSitesm2f3=seqSites;
    seqSitesm1delete1 = seqSites;
    ExprPredictor* predictor2 = new ExprPredictor(seqSitesb, bindingData, seqSites, seqLengths, exprData, 
        motifs, factorExprData, intFunc, coopMat, actIndicators, maxContact, repIndicators, repressionMat, 
        repressionDistThr, mmm, mmmr, indicator_bool3, anny, exprFile, seqSitesbot, seqSitesm1,seqSitesm2, 
        seqSitesf2 ,seqSitesbotf2, seqSitesm1f2 ,seqSitesm2f2, seqSitesf3, seqSitesbotf3,seqSitesm1f3,
         seqSitesm2f3 );  
    rval = readSequences( adamiFile, ExprPredictor::seqsy, ExprPredictor::seqNmes );
    // predictor2->objFuncborder(  par_init );
    // predictor2->train( predictor->getPar(), rng);
    predictor->objFuncborder(  par_init );
    predictor->train( predictor->getPar(), rng);
    // predictor2->train4( predictor2->getPar() );
    ofstream to("oData/ot.txt");
    // try {
    //     predictor->printFile3c(to,predictor2->getPar(),*predictor);
    // } catch (...) {
    //     predictor->printPar(predictor2->getPar());
    // } 
    try {
        predictor->printFile3c(to,predictor->getPar(),*predictor);
    } catch (...) {
        predictor->printPar(predictor->getPar());
    }
    to.close();
    // exit(0);
    ofstream fout( outFile.c_str() );
    if ( !fout ) {
        cerr << "Cannot open file " << outFile << endl;
        exit( 1 );
    }
    fout << "Rows\t" << condNames << endl;
    for ( int i = 0; i < nSeqs; i++ ) {
        vector< double > observedExprs = exprData2.getRow( i );
        vector< double > targetExprsm1;
        vector< double > dorsalExprs = factorExprData.getRow( 2 );
        predictor->predict(seqSitesm1[ i ], seqLengths[i], targetExprsm1 );
        fout << seqNames[i] << '\t';
        fout << observedExprs << endl;      
        fout << seqNames[i];
        for ( int j = 0; j < nConds; j++ ) {
            if(j==0){ fout << "\t" <<  targetExprsm1[j]+.0001; continue;} 
        fout << "\t" <<  targetExprsm1[j];       
    }
    fout << endl;
    fout << "Dl" << '\t';
    fout << dorsalExprs << endl;      
}
    fout.close();
return 0;	
}
