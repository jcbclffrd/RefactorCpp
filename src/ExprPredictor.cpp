#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include "ExprPredictor.h"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>
#include <sstream>
bool isNt( int a )
{
    if ( a < 0 || a > 3 ) return false;
    else return true;
}
int complement( int a )
{
    assert( a >= 0 && a < ALPHABET_SIZE );
    if ( a == 0 ) return 3;
    if ( a == 1 ) return 2;
    if ( a == 2 ) return 1;
    if ( a == 3 ) return 0;
    if ( a == MISSING ) return MISSING;
    if ( a == GAP ) return GAP;
}
int symbolToInt( char c )
{
    char upper = tolower( c );
    for ( int i = 0; i < ALPHABET_SIZE; i++ ) {
        if ( ALPHABET[ i ] == upper ) return i;
    }
    return -1;
}
char strand2char( bool strand )
{
    if ( strand ) return '+';
    else return '-';
}
bool char2strand( char c )
{
    assert( c == '+' || c == '-' );
    if ( c == '+' ) return true;
    else return false;
}
vector< double > createNtDistr( double gcContent )
{
    assert( gcContent >= 0 && gcContent <= 1.0 );
    vector< double > freqs( 4 );
    freqs[0] = ( 1.0 - gcContent ) / 2.0;
    freqs[1] = gcContent / 2.0;
    freqs[2] = freqs[1];
    freqs[3] = freqs[0];
    return freqs;
}
Sequence::Sequence( const string& str )
{
    for ( int i = 0; i < str.size(); i++ ) {
        int nt = symbolToInt( str[ i ] );	// could be a NNN or gap
        if ( nt >= 0 && nt < ALPHABET_SIZE ) {
            nts.push_back( nt );
        } else {
            cerr << "Illegal symbol: " << nt << " in " << str << endl;
            exit( 0 );
        }
    }
}
Sequence::Sequence( const Sequence& other, int start, int length, bool strand )
{
    assert( start >= 0 && length >= 0 && ( start + length ) <= other.size() );
    for ( int i = 0; i < length; i++ ) {
        if ( strand ) {	nts.push_back( other[ start + i ] ); }
        else { nts.push_back( complement( other[ start + length - 1 - i ] ) ); }
    }
}
int Sequence::push_back( int nt )
{
    assert( nt >= 0 && nt < ALPHABET_SIZE );
    nts.push_back( nt );
    return 0;
}
int Sequence::push_back( const Sequence& elem )
{
    for ( int i = 0; i < elem.size(); i++ ) push_back( elem[ i ] );
    return 0;
}
Sequence Sequence::compRevCompl() const
{
    return Sequence( *this, 0, size(), false );
}
void Sequence::getNtCounts( vector< int >& counts ) const
{
    counts.clear();
    for ( int i = 0; i < NBASES; i++ ) {
        counts.push_back( 0 );
    }
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] != GAP ) counts[ nts[ i ] ]++;
    }
}
bool Sequence::containsMissing() const
{
    for ( int i = 0; i < nts.size(); i++ ) {
        if ( nts[ i ] == MISSING ) return true;
    }
    return false;
}
int Sequence::load( const string& file, string& name, int format )
{
    vector< Sequence > seqs;
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );
    if ( rval == RET_ERROR ) return RET_ERROR;
    copy( seqs[ 0 ] );
    name = names[ 0 ];
    return rval;
}
int Sequence::load( const string& file, int format )
{
    string name;
    int rval = load( file, name, format );
    return rval;
}
ostream& operator<<( ostream& os, const Sequence& seq )
{
    for ( int i = 0; i < seq.size(); i++ ) {
        os << ALPHABET[ seq[ i ] ];
    }
    return os;
}
int readSequences( const string& file, vector< Sequence >& seqs, vector< string >& names, int format )
{
    if ( format != FASTA ) { return RET_ERROR; }
    seqs.clear();
    names.clear();
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }
    string line;
    Sequence seq;
    if ( format == FASTA ) {
        while ( getline( fin, line ) ) {
            if ( line[ 0 ] == '>' ) {
                if ( seq.size() ) {
                    seqs.push_back( seq );
                    seq.clear();
                }
                stringstream ss( line.substr( 1 ) );
                string name;
                ss >> name;
                names.push_back( name );
            } else {
                int start = line.find_first_not_of( " \t\r" );
                int last = line.find_last_not_of( " \t\r" );
                if ( start == string::npos || last == string::npos ) continue;
                for ( int i = start; i <= last; i++ ) {
                    int nt = symbolToInt( line[ i ] );	// could be a NNN or gap
                    if ( nt >= 0 && nt < ALPHABET_SIZE ) {
                        seq.push_back( nt );
                    } else {
                        return RET_ERROR;
                    }
                }
            }
        }
        if( seq.size() ) seqs.push_back( seq );
        return 0;
    }
}
int readSequences( const string& file, vector< Sequence >& seqs, int format )
{
    vector< string > names;
    int rval = readSequences( file, seqs, names, format );
    return rval;
}
int writeSequences( const string& file, const vector< Sequence >& seqs, const vector< string >& names, int format )
{
    assert( seqs.size() == names.size() );
    if ( format != FASTA ) { return RET_ERROR; }
    ofstream fout( file.c_str() );
    if ( format == FASTA ) {
        for ( int i = 0; i < seqs.size(); i++ ) {
            fout << ">" << names[ i ] << endl;
            fout << seqs[ i ] << endl;
        }
    }
    return 0;
}
int writeSequences( const string& file, const vector< Sequence >& seqs, int format )
{
    vector< string > names;
    for ( int i = 0; i < seqs.size(); i++ ) {
        char buffer[ 10 ];
        sprintf( buffer, "%i", i );
        names.push_back( string( buffer ) );
    }
    return writeSequences( file, seqs, names, format );
}
Matrix compWtmx( const Matrix& countMatrix, double pseudoCount )
{
    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
    int l = countMatrix.nRows();		// l: the length of motif
    Matrix pwm( l, 4 );
    for ( int i = 0; i < l; i++ ) {
        double n = 0;       // total counts at this position
        for ( int j = 0; j < 4; j++ ) {
            n += countMatrix( i, j );
        }
        for ( int j = 0; j < 4; j++ ) {
            pwm( i, j ) = ( countMatrix( i, j ) + pseudoCount ) / ( n + 4.0 * pseudoCount );
        }
    }
    return pwm;
}
Motif::Motif( const Matrix& _pwm, const vector< double >& _background ) : pwm( _pwm ), background( _background ), LLRMat( pwm.nRows(), 4 )
{
    assert( background.size() == 4 );
    init();
}
Motif::Motif( const Matrix& countMatrix, double pseudoCount, const vector< double >& _background ) : background( _background ), LLRMat( countMatrix.nRows(), 4 )
{
    assert( background.size() == 4 );
    pwm = compWtmx( countMatrix, pseudoCount );
    init();
}
double Motif::LLR( const Sequence& elem ) const
{
    int l = pwm.nRows();
    if ( elem.size() != l ) return GSL_NEGINF;
    if ( elem.containsMissing() ) return GSL_NEGINF;
    double result = 0;
    for ( int i = 0; i < l; i++ ) {
        result += LLRMat( i, elem[ i ] ); // LLR matrix, M(i,b) = log( f_i(b) / p(b) ), where f_i(b) is the frequency of b at position i of PWM, and p(b) the frequency of b at the background
    }
    return result;
}
double Motif::energy( const Sequence& elem ) const
{
    return ( -LLR( elem ) +   maxLLR );	  /// this is delta G,  Gmaxsite - Gsite, expratio = exp( - dG ) = exp ( LLR- LLRmax), since LLRmax =10, and LLR < 10 it follows
}
void Motif::sample( const gsl_rng* rng, Sequence& elem, bool strand ) const
{
    assert( rng != NULL );
    int l = pwm.nRows();
    Sequence sampleElem;
    for ( int i = 0; i < l; i++ ) {
        vector< double > distr = pwm.getRow( i );
        int nt = sampleMul( rng, distr );
        sampleElem.push_back( nt );
    }
    if ( strand == 0 ) elem = sampleElem.compRevCompl();
    else elem = sampleElem;
}
int Motif::load( const string& file, const vector< double >& background, string& name )
{
    vector< Motif > motifs;
    vector< string > names;
    int rval = readMotifs( file, background, motifs, names );
    if ( rval == RET_ERROR ) return RET_ERROR;
    copy( motifs[ 0 ] );
    name = names[ 0 ];
    return rval;
}
int Motif::load( const string& file, const vector< double >& background )
{
    string name;
    int rval = load( file, background, name );
    return rval;
}
ostream& operator<<( ostream& os, const Motif& motif )
{
    os << motif.pwm;
    return os;
}
void Motif::init()
{
    int l = pwm.nRows();
    for ( int i = 0; i < l; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
            LLRMat( i, j ) = log( pwm( i, j ) / background[ j ] );
        }
    }
    for ( int i = 0; i < l; i++ ) {
        int b_max;
        max( pwm.getRow( i ), b_max );
        maxSite.push_back( b_max );
    }
    maxLLR = 0;
    for ( int i = 0; i < l; i++ ) {
        maxLLR += LLRMat( i, maxSite[ i ] );
    }
}
Matrix countmatrixS( const string& file )
{
    ifstream seq_file( file.c_str() );
    if ( !seq_file ) {
    }
    string temp;
    vector <string> seq_name;
    vector <string> seq;
    string templ;
    int tem = 0;
    int tempmax = 0;
    seq.clear();
    while(!seq_file.eof()){
        temp = "";
        getline(seq_file, temp);
        if(temp.length() == 0){
            break;
        }
        string name (temp, 1, temp.length() - 1);
        seq_name.push_back(name);
        getline(seq_file, temp);
        seq.push_back(temp);
        tem = temp.length();
        if( tem > tempmax) tempmax = tem;
    }
    int length = tempmax ;  //maxstring(file.c_str());
    Matrix m(4,length,0);
    for( int j = 0; j < seq.size(); j++ ){
        Sequence readseq(seq[j]);
        for( int i = 0;  i< readseq.size(); i++)  // this for loop is miss-counting by one? maybe overloaded size method?
        {
            if(readseq[i] == 0)            //const char ALPHABET[] = { 'a', 'c', 'g', 't', 'N', '-' }
            { m(0,i) = m(0,i) + 1;}
            if(readseq[i] == 1)
            { m(1,i) = m(1,i) + 1;}
            if(readseq[i] == 2)
            { m(2,i) = m(2,i) + 1;}
            if(readseq[i] == 3)
            { m(3,i) = m(3,i) + 1;}
        }
    } // for j
    cout << " m "  << endl << m << endl;
    Matrix countMatrix = m.transpose();
    cout << "countmat" << countMatrix << endl;
    double pseudoCount =.1;
    assert( countMatrix.nCols() == 4 && pseudoCount >= 0 );
    int l = countMatrix.nRows();		// l: the length of motif
    Matrix   pwm = compWtmx( countMatrix, pseudoCount ) ;
    return pwm;
}
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs, vector< string >& names )
{
    ifstream fin( file.c_str() );
    if ( !fin ) { cerr << "Cannot open" << file << endl; exit( 1 ); }
    motifs.clear();
    names.clear();
    string line;
    do {
        getline( fin, line );
        if ( line[ 0 ] != '>' ) continue;
        int MAX_SIZE = 100;
        char lineStr[ MAX_SIZE ];
        strcpy( lineStr, ( line.substr( 1 ) ).c_str() );
        char *name, *lengthStr, *pseudoCountStr;
        name = strtok( lineStr, " \t" );
        lengthStr = strtok( NULL, " \t" );
        pseudoCountStr = strtok( NULL, " \t" );
        int length;
        double pseudoCount;
        if ( lengthStr ) length = atoi( lengthStr );
        else { return RET_ERROR; }
        if ( pseudoCountStr ) pseudoCount = atof( pseudoCountStr );
        else pseudoCount = PSEUDO_COUNT;
        Matrix countMat( length, NBASES );
        for ( int i = 0; i < length; ++i ) {
            for ( int j = 0; j < NBASES; ++j ) {
                fin >> countMat( i, j );
            }
        }
        names.push_back( string( name ) );
        motifs.push_back( Motif( countMat, pseudoCount, background ) );
    } while ( !fin.eof() );
    return 0;
}
int readMotifs( const string& file, const vector< double >& background, vector< Motif >& motifs )
{
    vector< string > names;
    return readMotifs( file, background, motifs, names );
}
ostream& operator<<( ostream& os, const Site& site )
{
    char strandChar = site.strand ? '+' : '-';
    os <<"\n"<< site.start + 1 << "\t" << strandChar << "\t" << site.factorIdx << "\t" << site.energy << "\t" << site.wtRatio;
    return os;
}
bool siteOverlap( const Site& a, const Site& b, const vector< Motif >& motifs )
{
    if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
    if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;
    return true;
}
bool siteOverlap2( const Site& a, const Site& b, const vector< Motif >& motifs )
{
    if ( a.start + motifs[ a.factorIdx ].length() <= b.start ) return false;
    if ( b.start + motifs[ b.factorIdx ].length() <= a.start ) return false;
    if( a.start == b.start ) return false;
    return true;
}
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, vector< string >& names, bool readEnergy )
{
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    sites.clear();
    names.clear();
    SiteVec currVec;
    int nrecords = 0;       // number of ">" read so far
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        if ( line.substr( 0, 1 ) == ">" ) {
            stringstream ss( line.substr( 1 ) );
            string name;
            ss >> name;
            names.push_back( name );
            nrecords++;
            if ( nrecords > 1 ) {
                sites.push_back( currVec );
                currVec.clear();
            }
        } else {
            int start;
            char strandChar;
            string factor;
            double energy = 0;
            stringstream ss( line );
            ss >> start >> strandChar >> factor;
            if ( readEnergy ) ss >> energy;
            bool strand = strandChar == '+' ? 1 : 0;
            map<string, int>::const_iterator iter = factorIdxMap.find( factor );
            currVec.push_back( Site( start - 1, strand, iter->second , energy ) );
        }
    }
    sites.push_back( currVec );
    return 0;
}
int readSites( const string& file, const map< string, int >& factorIdxMap, vector< SiteVec >& sites, bool readEnergy )
{
    vector< string > names;
    return readSites( file, factorIdxMap, sites, names, readEnergy );
}
int SeqAnnotator::annot( const Sequence& seq, SiteVec& sites ) const
{
    sites.clear();
    for ( int i = 0; i < seq.size(); i++ ) {
        for ( int k = 0; k < motifs.size(); k++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) continue;
            double energy;
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrs[ k ] ) {  //1115 changed from <= to >=
            sites.push_back( Site( i, 1, k, energy ) );
        }
        Sequence rcElem( seq, i, l, 0 );
        energy = motifs[ k ].energy( rcElem );
        if ( energy <= energyThrs[ k ] ) {//1115 changed from <= to >= , 1129 changed back
        sites.push_back( Site( i, 0, k, energy ) );
    }
}
}
return sites.size();
}
bool siteSortPredicate(const Site& d1, const Site& d2)
{
    return d1.start < d2.start;
}
int SeqAnnotator::annoty4( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names)
{
    vector< double > cons = f.getCol( column );
    tsites.clear();
    Site t;
    vector< double > ff(3);  // try this without initializing
    typedef  vector< Site > sitestt;
    sitestt te;
    vector< sitestt > sitest( motifs.size() );    // sitestemp
    vector< sitestt > sitesp( motifs.size() );    // sitespermanent
    for ( int k = 0; k < 2; k++ ) {
        for ( int i = 0; i < seq.size(); i++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) break ;
            double energy;
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 1, k, energy ) );
            }
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 0, k, energy ) );
            }
        }
        if( k == 0 ) {
            double mine =100;
            int gbest;
            int counter=1;
            sitestt  s = sitest[0];
            tsites.clear();
            while( counter > 0 && s.size() > 0 ){
                counter=0;
                Site tempsite;
                for ( int i = 0; i < s.size(); i++ ) {
                    if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
                    mine = s[i].energy;
                    tempsite = s[i];
                    gbest = i;
                    counter++;
                }
            }
            tsites.push_back (tempsite);
            sitestoverlap( s, tsites );
        } // while
        std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
        p.predictOcc( tsites , 5, cons, ff);
        if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
        return 1;
    }
    else{ tsites.clear();}
}  // if k
}  // for k
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
    ktempindex = k;
    for ( int i = 0; i < sitest[1].size(); i++ ) {
        sitesorderedtemp.push_back( sitest[0][k] );
        sitesorderedtemp.push_back( sitest[1][i] );
        std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
        Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
        sitesorderedtemp.clear();
        itempindex=i;
        if( Zbest < Ztemp ) {
            Zbest = Ztemp;
            ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
            kbest = ktempindex;
        }// if
    }// for i
}
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] );
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] );
vector< double  > sitest_tw_Z( 0 );
int N=1;
double previous = 0;
double current = 0;
do {
    tsites.clear();
    for ( int k = 0; k < 2; k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
    p.predictOcc( tsites, 5, cons, ff);
    previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
    if( previous > .5 )
    { break ; }
    N++;
    tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal
    for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
    sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
    std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
    sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
    sitesorderedtemp.clear();
}   // for i
int ind = maxZ( sitest_tw_Z );
double maxZscore = maxZs( sitest_tw_Z );
sitesp[1].push_back( sitest[1][ ind ] );
sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;
tsites.clear();
for ( int k = 0; k < 2; k++ ) {
    for( int i = 0; i < sitesp[k].size(); i++){
        tsites.push_back( sitesp[k][i] );
    }
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];
if( current > .5 ) break;
} while ( abs(previous - current ) > .25 );
return tsites.size();
}
int SeqAnnotator::annotydorsalold(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names)
{
    vector< double > cons = f.getCol( column );
    tsites.clear();
    Site t;
    vector< double > ff(3);  // try this without initializing
    if( cons[2] ==0 ) {
        typedef  vector< Site > sitestt;
        sitestt te;
        vector< sitestt > sitest( motifs.size() );    // sitestemp
        vector< sitestt > sitesp( motifs.size() );    // sitespermanent
        for ( int k = 0; k < 2; k++ ) {
            for ( int i = 0; i < seq.size(); i++ ) {
                int l = motifs[ k ].length();
                if ( i + l > seq.size() ) break ;
                double energy;
                Sequence elem( seq, i, l, 1 );
                energy = motifs[ k ].energy( elem );
                if ( energy <= energyThrs[ k ] ) {
                    sitest[k].push_back( Site( i, 1, k, energy ) );
                }
                Sequence rcElem( seq, i, l, 0 );
                energy = motifs[ k ].energy( rcElem );
                if ( energy <= energyThrs[ k ] ) {
                    sitest[k].push_back( Site( i, 0, k, energy ) );
                }
            }
        }  // for k
        SiteVec sitesorderedtemp(0);
        double Ztemp=0;
        int ktempindex=0;
        int itempindex=0;
        double Zbest=0;
        int ibest =0;
        int  kbest=0;
        for ( int k = 0; k < sitest[0].size(); k++ ) {
            ktempindex = k;
            for ( int i = 0; i < sitest[1].size(); i++ ) {
                sitesorderedtemp.push_back( sitest[0][k] );
                sitesorderedtemp.push_back( sitest[1][i] );
                std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
                Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
                sitesorderedtemp.clear();
                itempindex=i;
                if( Zbest < Ztemp ) {
                    Zbest = Ztemp;
                    ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
                    kbest = ktempindex;
                }// if
            }// for i
        }
        sitesp[0].push_back( sitest[0][kbest] );
        sitestoverlap( sitest[0], sitesp[0] );
        sitesp[1].push_back( sitest[1][ibest] );
        sitestoverlap( sitest[1], sitesp[1] );
        vector< double  > sitest_tw_Z( 0 );
        int N=1;
        double previous = 0;
        double current = 0;
        do {
            tsites.clear();
            for ( int k = 0; k < 2; k++ ) {
                for( int i = 0; i < sitesp[k].size(); i++){
                    tsites.push_back( sitesp[k][i] );
                }
            }
            std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
            p.predictOcc( tsites, 5, cons, ff);
            previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
            if( previous > .5 )
            { break ; }
            N++;
            tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal
            for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
            sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
            sitesorderedtemp.clear();
        }   // for i
        int ind = maxZ( sitest_tw_Z );
        double maxZscore = maxZs( sitest_tw_Z );
        sitesp[1].push_back( sitest[1][ ind ] );
        sitestoverlap( sitest[1], sitesp[1] );
        if (N > 6) break;
        tsites.clear();
        for ( int k = 0; k < 2; k++ ) {
            for( int i = 0; i < sitesp[k].size(); i++){
                tsites.push_back( sitesp[k][i] );
            }
        }
        std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
        p.predictOcc( tsites, 5, cons, ff);
        current = ff[0];
        if( current > .5 ) break;
    } while ( abs(previous - current ) > .25 );
    return tsites.size();
} //if c=0
else {
    typedef  vector< Site > sitestt;
    vector< sitestt > sitest( motifs.size() );    // sitestemp
    vector< sitestt > sitesp( motifs.size() );    // sitespermanent
    for ( int k = 0; k < motifs.size(); k++ ) {
        for ( int i = 0; i < seq.size(); i++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) break ;
            double energy;
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 1, k, energy ) );
            }
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 0, k, energy ) );
            }
        }
        if( k == 0 ) {
            double mine =100;
            int gbest;
            int counter=1;
            sitestt  s = sitest[0];
            tsites.clear();
            while( counter > 0 && s.size() > 0 ){
                counter=0;
                Site tempsite;
                for ( int i = 0; i < s.size(); i++ ) {
                    if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
                    mine = s[i].energy;
                    tempsite = s[i];
                    gbest = i;
                    counter++;
                }
            }
            if( counter != 0) {
                tsites.push_back (tempsite);
                sitestoverlap( s, tsites );
            }
        } // while
        std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
        p.predictOcc( tsites , 5, cons, ff);
        if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
        cout << tsites << endl;
        s.clear();
        return 1;
    }
    else{ s.clear(); tsites.clear();}
}  // if k
} //for k
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
    ktempindex = k;
    for ( int i = 0; i < sitest[1].size(); i++ ) {
        sitesorderedtemp.push_back( sitest[0][k] );
        sitesorderedtemp.push_back( sitest[1][i] );
        std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
        Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
        sitesorderedtemp.clear();
        itempindex=i;
        if( Zbest < Ztemp ) {
            Zbest = Ztemp;
            ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
            kbest = ktempindex;
        }// if
    }// for i
}
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
    ktempindex = k;
    for ( int i = k+1; i < sitest[0].size(); i++ ) {
        sitesorderedtemp.push_back( sitest[0][k] );
        sitesorderedtemp.push_back( sitest[1][i] );
        std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
        Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
        sitesorderedtemp.clear();
        itempindex=i;
        if( Zbest < Ztemp ) {
            Zbestdd = Ztemp;
            ibestdd = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
            kbestdd = ktempindex;
        }// if
    }// for i
}
if( Zbestdd < Zbest ) {
    sitesp[0].push_back( sitest[0][kbest] );
    sitestoverlap( sitest[0], sitesp[0] );
    p.predictOcc( sitesp[0] , 5, cons, ff);
    if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
    tsites= sitesp[0];
    return 1; }
    sitesp[1].push_back( sitest[1][ibest] );
    sitestoverlap( sitest[1], sitesp[1] );
}
else{                                     // if 2 dorsal score higher than a dt pair then push dorsals
sitesp[0].push_back( sitest[0][kbestdd] );
sitestoverlap( sitest[0], sitesp[0] );
sitesp[0].push_back( sitest[0][ibestdd] );
sitestoverlap( sitest[0], sitesp[0] );
p.predictOcc( sitesp[0] , 5, cons, ff);
if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
tsites= sitesp[0];
return 1; }
}
sitesp[2].push_back (siteMax( sitest[2] ));
sitestoverlap( sitest[2], sitesp[2] );
int N = 1; // this is the number of sites
double previous;
double current;
do {
    vector< double > ff(3); // focc passed by ref to predictocc.
    tsites.clear();
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    N++;
    vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
    SiteVec sitesorderedtemp = tsites;
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitest[k].size(); i++){
            sitesorderedtemp = tsites;
            sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
            sitesorderedtemp.clear();
        }   // for i
        vector< double > uniq = sitest_allk_alli_Z[k];  // maybe multiply this by 1000 and use int for unique (thereby getting the right truncation error.
        std::unique(uniq.begin(),uniq.end() );
        if ( uniq.size() == 1 ) {
            continue;
        }
        int ind = maxZ( sitest_allk_alli_Z[k] );
        if ( cons[k] != 0 ) {
            sitesp[k].push_back( sitest[k][ ind ] );
            sitestoverlap( sitest[k], sitesp[k] );
        }
    }  // for k
    tsites.clear();
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
    p.predictOcc( tsites, 5, cons, ff);
    current = ff[0];  // using foccdl in the while condition is out of scope.  why?
    if( current > .5 )
    { break ; }
    if(N ==6) {
        break;
    }
} while( N < 3 ) ;//(abs(previous - current) > .25 ) );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
    for( int i = 0; i < sitesp[k].size(); i++){
        tsites.push_back( sitesp[k][i] );
    }
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
return sitest.size();
}// else
}
int SeqAnnotator::annotydorsal(const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names)
{
    vector< double > cons = f.getCol( column );
    tsites.clear();
    Site t;
    vector< double > ff(3);  // try this without initializing
    if( cons[2] ==0 ) {
        typedef  vector< Site > sitestt;
        sitestt te;
        vector< sitestt > sitest( motifs.size() );    // sitestemp
        vector< sitestt > sitesp( motifs.size() );    // sitespermanent
        for ( int k = 0; k < 2; k++ ) {
            for ( int i = 0; i < seq.size(); i++ ) {
                int l = motifs[ k ].length();
                if ( i + l > seq.size() ) break ;
                double energy;
                Sequence elem( seq, i, l, 1 );
                energy = motifs[ k ].energy( elem );
                if ( energy <= energyThrs[ k ] ) {
                    sitest[k].push_back( Site( i, 1, k, energy ) );
                }
                Sequence rcElem( seq, i, l, 0 );
                energy = motifs[ k ].energy( rcElem );
                if ( energy <= energyThrs[ k ] ) {
                    sitest[k].push_back( Site( i, 0, k, energy ) );
                }
            } // for i
            if( k == 0 ) {
                for ( int i = 0; i < seq.size(); i++ ) {
                    int l = d2.length();
                    if ( i + l > seq.size() ) break ;
                    double energy;
                    Sequence elem( seq, i, l, 1 );
                    energy = d2.energy( elem );
                    if ( energy <= energyThrs[ 0 ] ) {
                        sitest[0].push_back( Site( i, 1, 0, energy ) );
                    }
                    Sequence rcElem( seq, i, l, 0 );
                    energy = d2.energy( rcElem );
                    if ( energy <= energyThrs[ 0] ) {
                        sitest[0].push_back( Site( i, 0, 0, energy ) );
                    }
                }  // for i
                double mine =100;
                int gbest;
                int counter=1;
                sitestt  s = sitest[k];
                tsites.clear();
                while( counter > 0 && s.size() > 0 ){
                    counter=0;
                    Site tempsite;
                    for ( int i = 0; i < s.size(); i++ ) {
                        if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
                        mine = s[i].energy;
                        tempsite = s[i];
                        gbest = i;
                        counter++;
                    }
                }
                if( counter != 0) {
                    tsites.push_back (tempsite);
                    sitestoverlap( s, tsites );
                }
            } // while
            mine =100;
            gbest = 100;
            counter=1;
            while( counter > 0 && s.size() > 0 ){  // in case d2 also has a set of motifs repeat while loop, making sure not to double count first motif of dl.
            counter=0;
            Site tempsite;
            for ( int i = 0; i < s.size(); i++ ) {
                if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
                mine = s[i].energy;
                tempsite = s[i];
                gbest = i;
                counter++;
            }
        }
        if( counter != 0) {
            tsites.push_back (tempsite);
            sitestoverlap( s, tsites );
        }
    } // while
    s.clear();
} // if k == 0
else{
    double mine =100;
    int gbest;
    int counter=1;
    sitestt  s = sitest[k];
    while( counter > 0 && s.size() > 0 ){
        counter=0;
        Site tempsite;
        for ( int i = 0; i < s.size(); i++ ) {
            if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
            mine = s[i].energy;
            tempsite = s[i];
            gbest = i;
            counter++;
        }
    }
    if( counter != 0) {
        tsites.push_back (tempsite);
        sitestoverlap( s, tsites );
    }
} // while
s.clear();
} // else
} //for k
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites , 5, cons, ff);
if ( ff[0] >.5) {
    return 1;
}
else{ tsites.clear();}
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
    ktempindex = k;
    for ( int i = 0; i < sitest[1].size(); i++ ) {
        sitesorderedtemp.push_back( sitest[0][k] );
        sitesorderedtemp.push_back( sitest[1][i] );
        std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
        Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
        sitesorderedtemp.clear();
        itempindex=i;
        if( Zbest < Ztemp ) {
            Zbest = Ztemp;
            ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
            kbest = ktempindex;
        }// if
    }// for i
} // for k
sitesp[0].push_back( sitest[0][kbest] );
sitestoverlap( sitest[0], sitesp[0] );
sitesp[1].push_back( sitest[1][ibest] );
sitestoverlap( sitest[1], sitesp[1] );
vector< double  > sitest_tw_Z( 0 );
int N=1;
double previous = 0;
double current = 0;
do {
    tsites.clear();
    for ( int k = 0; k < 2; k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
    p.predictOcc( tsites, 5, cons, ff);
    previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
    if( previous > .5 )
    { break ; }
    N++;
    tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal; otherwise place for loop over all dorsal sites here( changing alg from n+n to n*n, where n is the number of sites
    for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
    sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
    std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
    sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
    sitesorderedtemp.clear();
}   // for i
int ind = maxZ( sitest_tw_Z );
double maxZscore = maxZs( sitest_tw_Z );
sitesp[1].push_back( sitest[1][ ind ] );
sitestoverlap( sitest[1], sitesp[1] );
if (N > 6) break;
tsites.clear();
for ( int k = 0; k < 2; k++ ) {
    for( int i = 0; i < sitesp[k].size(); i++){
        tsites.push_back( sitesp[k][i] );
    }
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites, 5, cons, ff);
current = ff[0];
if( current > .5 ) break;
}
while ( abs(previous - current ) > .25 );
return tsites.size();
} //if( cons[2] ==0 )
else {
    typedef  vector< Site > sitestt;
    vector< sitestt > sitest( motifs.size() );    // sitestemp
    vector< sitestt > sitesp( motifs.size() );    // sitespermanent
    for ( int k = 0; k < motifs.size(); k++ ) {
        for ( int i = 0; i < seq.size(); i++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) break ;
            double energy;
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 1, k, energy ) );
            }
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 0, k, energy ) );
            }
        }
        if( k == 0 ) {
            for ( int i = 0; i < seq.size(); i++ ) {
                int l = d2.length();
                if ( i + l > seq.size() ) break ;
                double energy;
                Sequence elem( seq, i, l, 1 );
                energy = d2.energy( elem );
                if ( energy <= energyThrs[ 0 ] ) {
                    sitest[0].push_back( Site( i, 1, 0, energy ) );
                }
                Sequence rcElem( seq, i, l, 0 );
                energy = d2.energy( rcElem );
                if ( energy <= energyThrs[ 0] ) {
                    sitest[0].push_back( Site( i, 0, 0, energy ) );
                }
            }
            double mine =100;
            int gbest;
            int counter=1;
            int count=0;
            sitestt  s = sitest[k];
            tsites.clear();
            while( counter > 0 && s.size() > 0 ){
                counter=0;
                Site tempsite;
                for ( int i = 0; i < s.size(); i++ ) {
                    if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
                    mine = s[i].energy;
                    tempsite = s[i];
                    gbest = i;
                    counter++;
                    count++;
                }
            }
            if( counter != 0) {
                tsites.push_back (tempsite);
                sitestoverlap( s, tsites );
            }
        } // while
        s.clear();
    }  // if k
    else{
        if( sitest[k].size() == 0 ) {continue; }
        double mine =100;
        int gbest;
        int counter=1;
        sitestt  s = sitest[k];
        while( counter > 0 && s.size() > 0 ){
            counter=0;
            Site tempsite;
            for ( int i = 0; i < s.size(); i++ ) {
                if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
                mine = s[i].energy;
                tempsite = s[i];
                gbest = i;
                counter++;
            }
        }
        if( counter != 0) {
            tsites.push_back (tempsite);
            sitestoverlap( s, tsites );
        }
    } // while
    s.clear();
}  // else
} //for k
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
p.predictOcc( tsites , 5, cons, ff);
if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
return 1;
}
else{ tsites.clear();}
SiteVec sitesorderedtemp(0);
double Ztemp=0;
int ktempindex=0;
int itempindex=0;
double Zbest=0;
int ibest =0;
int  kbest=0;
for ( int k = 0; k < sitest[0].size(); k++ ) {
    ktempindex = k;
    for ( int i = 0; i < sitest[1].size(); i++ ) {
        sitesorderedtemp.push_back( sitest[0][k] );
        sitesorderedtemp.push_back( sitest[1][i] );
        std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
        Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
        sitesorderedtemp.clear();
        itempindex=i;
        if( Zbest < Ztemp ) {
            Zbest = Ztemp;
            ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
            kbest = ktempindex;
        }// if
    }// for i
} // for k
double Zbestdd=0;
int ibestdd =0;
int kbestdd = 0;
if( sitest[0].size() > 1 ) {
    for ( int k = 0; k < sitest[0].size(); k++ ) {
        ktempindex = k;
        for ( int i = k+1; i < sitest[0].size(); i++ ) {
            sitesorderedtemp.push_back( sitest[0][k] );
            sitesorderedtemp.push_back( sitest[1][i] );
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
            sitesorderedtemp.clear();
            itempindex=i;
            if( Zbest < Ztemp ) {
                Zbestdd = Ztemp;
                ibestdd = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
                kbestdd = ktempindex;
            }// if
        }// for i
    } //for k
    if( Zbestdd < Zbest ) {
        sitesp[0].push_back( sitest[0][kbest] );
        sitestoverlap( sitest[0], sitesp[0] );
        p.predictOcc( sitesp[0] , 5, cons, ff);
        if ( ff[0] >.5) { // cout << " ff[0] " << ff[0] << endl;
        tsites= sitesp[0];
        return 1; }
    } // if zbes < zbe
    else{                                     // if 2 dorsal score higher than a dt pair then push dorsals
    sitesp[0].push_back( sitest[0][kbestdd] );
    sitestoverlap( sitest[0], sitesp[0] );
    sitesp[0].push_back( sitest[0][ibestdd] );
    sitestoverlap( sitest[0], sitesp[0] );
    p.predictOcc( sitesp[0] , 5, cons, ff);
    if ( ff[0] >.5) { // cout << " ff[0] " << ff[0] << endl;
    tsites= sitesp[0];
    return 1; }
}
}// if sitest0 > 1
int N = 1; // this is the number of sites
double previous;
double current;
do {
    vector< double > ff(3); // focc passed by ref to predictocc.
    tsites.clear();
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    N++;
    vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
    SiteVec sitesorderedtemp = tsites;
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitest[k].size(); i++){
            sitesorderedtemp = tsites;
            sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
            sitesorderedtemp.clear();
        }   // for i
        if(sitest[k].size() > 0  ) {
            int ind = maxZ( sitest_allk_alli_Z[k] );
            if ( cons[k] != 0 ) {
                sitesp[k].push_back( sitest[k][ ind ] );
                sitestoverlap( sitest[k], sitesp[k] );
            }
        }
    }  // for k
    tsites.clear();
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
    vector< double > ff2(3);
    p.predictOcc( tsites, 5, cons, ff2);
    current = ff[0];  // using foccdl in the while condition is out of scope.  why?
    if( current > .5 )
    { break ; }
    if(N ==6) {
        break;
    }
} while( N < 3 ) ;//(abs(previous - current) > .25 ) );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
    for( int i = 0; i < sitesp[k].size(); i++){
        tsites.push_back( sitesp[k][i] );
    }
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
return sitest.size();
}// else
}
int SeqAnnotator::annoty2( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names)
{
    vector< double > cons = f.getCol( column );
    tsites.clear();
    Site t;
    vector< double > ff(3);  // try this without initializing
    cout << "neuro" << endl;
    if(cons[2] ==0 ) {
        typedef  vector< Site > sitestt;
        sitestt te;
        vector< sitestt > sitest( motifs.size() );    // sitestemp
        vector< sitestt > sitesp( motifs.size() );    // sitespermanent
        for ( int k = 0; k < 2; k++ ) {
            for ( int i = 0; i < seq.size(); i++ ) {
                int l = motifs[ k ].length();
                if ( i + l > seq.size() ) break ;
                double energy;
                Sequence elem( seq, i, l, 1 );
                energy = motifs[ k ].energy( elem );
                if ( energy <= energyThrs[ k ] ) {
                    sitest[k].push_back( Site( i, 1, k, energy ) );
                }
                Sequence rcElem( seq, i, l, 0 );
                energy = motifs[ k ].energy( rcElem );
                if ( energy <= energyThrs[ k ] ) {
                    sitest[k].push_back( Site( i, 0, k, energy ) );
                }
            }
            if( k == 0 ) {
                tsites.clear();
                tsites.push_back (siteMax( sitest[0] ));
                p.predictOcc( tsites , 5, cons, ff);
                if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
                return 1; }
                else{ tsites.clear();}
            }
        }
        SiteVec sitesorderedtemp(0);
        double Ztemp=0;
        int ktempindex=0;
        int itempindex=0;
        double Zbest=0;
        int ibest =0;
        int  kbest=0;
        for ( int k = 0; k < sitest[0].size(); k++ ) {
            ktempindex = k;
            for ( int i = 0; i < sitest[1].size(); i++ ) {
                sitesorderedtemp.push_back( sitest[0][k] );
                sitesorderedtemp.push_back( sitest[1][i] );
                std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
                Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
                sitesorderedtemp.clear();
                itempindex=i;
                if( Zbest < Ztemp ) {
                    Zbest = Ztemp;
                    ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
                    kbest = ktempindex;
                }// if
            }// for i
            cout << " not inside Zbest " << endl;
        }
        sitesp[0].push_back( sitest[0][kbest] );
        sitestoverlap( sitest[0], sitesp[0] );
        sitesp[1].push_back( sitest[1][ibest] );
        sitestoverlap( sitest[1], sitesp[1] );
        vector< double  > sitest_tw_Z( 0 );
        int N=1;
        double previous = 0;
        double current = 0;
        do {
            tsites.clear();
            for ( int k = 0; k < 2; k++ ) {
                for( int i = 0; i < sitesp[k].size(); i++){
                    tsites.push_back( sitesp[k][i] );
                }
            }
            std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
            p.predictOcc( tsites, 5, cons, ff);
            previous = ff[0];  // using foccdl in the while condition is out of scope.  why?
            if( previous > .5 )
            { break ; }
            N++;
            tsites.push_back( siteMax( sitest[0] ) );  // add another dorsal
            for( int i = 0; i < sitest[1].size(); i++){  // check all the rest of the twists sites
            sitesorderedtemp.push_back( sitest[1][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            sitest_tw_Z.push_back( p.predictZ( sitesorderedtemp,  cons) );
            sitesorderedtemp.clear();
        }   // for i
        int ind = maxZ( sitest_tw_Z );
        double maxZscore = maxZs( sitest_tw_Z );
        sitesp[1].push_back( sitest[1][ ind ] );
        sitestoverlap( sitest[1], sitesp[1] );
        if (N > 6) break;
        tsites.clear();
        for ( int k = 0; k < 2; k++ ) {
            for( int i = 0; i < sitesp[k].size(); i++){
                tsites.push_back( sitesp[k][i] );
            }
        }
        std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
        p.predictOcc( tsites, 5, cons, ff);
        current = ff[0];
        if( current > .5 ) break;
    } while ( abs(previous - current ) > .25 );
    return tsites.size();
}
else {
    typedef  vector< Site > sitestt;
    vector< sitestt > sitest( motifs.size() );    // sitestemp
    vector< sitestt > sitesp( motifs.size() );    // sitespermanent
    for ( int k = 0; k < motifs.size(); k++ ) {
        for ( int i = 0; i < seq.size(); i++ ) {
            int l = motifs[ k ].length();
            if ( i + l > seq.size() ) break ;
            double energy;
            Sequence elem( seq, i, l, 1 );
            energy = motifs[ k ].energy( elem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 1, k, energy ) );
            }
            Sequence rcElem( seq, i, l, 0 );
            energy = motifs[ k ].energy( rcElem );
            if ( energy <= energyThrs[ k ] ) {
                sitest[k].push_back( Site( i, 0, k, energy ) );
            }
        }
        if( k == 0 ) {
            tsites.clear();
            tsites.push_back (siteMax( sitest[0] ));
            p.predictOcc( tsites , 5, cons, ff);
            if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
            return 1; }
            else{ tsites.clear();}
        }  //if k
    } //for k
    SiteVec sitesorderedtemp(0);
    double Ztemp=0;
    int ktempindex=0;
    int itempindex=0;
    double Zbest=0;
    int ibest =0;
    int  kbest=0;
    for ( int k = 0; k < sitest[0].size(); k++ ) {
        ktempindex = k;
        for ( int i = 0; i < sitest[1].size(); i++ ) {
            sitesorderedtemp.push_back( sitest[0][k] );
            sitesorderedtemp.push_back( sitest[1][i] );
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
            sitesorderedtemp.clear();
            itempindex=i;
            if( Zbest < Ztemp ) {
                Zbest = Ztemp;
                ibest = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
                kbest = ktempindex;
            }// if
        }// for i
    }
    double Zbestdd=0;
    int ibestdd =0;
    int kbestdd = 0;
    for ( int k = 0; k < sitest[0].size(); k++ ) {
        ktempindex = k;
        for ( int i = k+1; i < sitest[0].size(); i++ ) {
            sitesorderedtemp.push_back( sitest[0][k] );
            sitesorderedtemp.push_back( sitest[1][i] );
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            Ztemp= p.predictZ( sitesorderedtemp,  cons) ;
            sitesorderedtemp.clear();
            itempindex=i;
            if( Zbest < Ztemp ) {
                Zbestdd = Ztemp;
                ibestdd = itempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
                kbestdd = ktempindex;
            }// if
        }// for i
    }
    if( Zbestdd < Zbest ) {
        sitesp[0].push_back( sitest[0][kbest] );
        sitestoverlap( sitest[0], sitesp[0] );
        p.predictOcc( sitesp[0] , 5, cons, ff);
        if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
        tsites= sitesp[0];
        return 1; }
        sitesp[1].push_back( sitest[1][ibest] );
        sitestoverlap( sitest[1], sitesp[1] );
    }
    else{                                     // if 2 dorsal score higher than a dt pair then push dorsals
    sitesp[0].push_back( sitest[0][kbestdd] );
    sitestoverlap( sitest[0], sitesp[0] );
    sitesp[0].push_back( sitest[0][ibestdd] );
    sitestoverlap( sitest[0], sitesp[0] );
    p.predictOcc( sitesp[0] , 5, cons, ff);
    if ( ff[0] >.5) {  //cout << " ff[0] " << ff[0] << endl;
    tsites= sitesp[0];
    return 1; }
}
int N = 1; // this is the number of sites
double previous;
double current;
do {
    vector< double > ff(3); // focc passed by ref to predictocc.
    tsites.clear();
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    N++;
    vector< vector < double > > sitest_allk_alli_Z( motifs.size() ); //this should be the same size as (sitest[k]), however the number of sites is : ( sitesp + sites[k] )
    SiteVec sitesorderedtemp = tsites;
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitest[k].size(); i++){
            sitesorderedtemp = tsites;
            sitesorderedtemp.push_back( sitest[k][i] )  ;  // reorder sitesorderedtemp, such that the pushed is ordered.
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            sitest_allk_alli_Z[k].push_back( p.predictZ( sitesorderedtemp,  cons) );
            sitesorderedtemp.clear();
        }   // for i
        vector< double > uniq = sitest_allk_alli_Z[k];  // maybe multiply this by 1000 and use int for unique (thereby getting the right truncation error.
        std::unique(uniq.begin(),uniq.end() );
        if ( uniq.size() == 1 ) {
            continue;
        }
        int ind = maxZ( sitest_allk_alli_Z[k] );
        if ( cons[k] != 0 ) {
            sitesp[k].push_back( sitest[k][ ind ] );
            sitestoverlap( sitest[k], sitesp[k] );
        }
    }  // for k
    tsites.clear();
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitesp[k].size(); i++){
            tsites.push_back( sitesp[k][i] );
        }
    }
    std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
    p.predictOcc( tsites, 5, cons, ff);
    current = ff[0];  // using foccdl in the while condition is out of scope.  why?
    if( current > .5 )
    { break ; }
    if(N ==6) {
        break;
    }
} while( N < 3 ) ;//(abs(previous - current) > .25 ) );  // no need to say ff[0] < 1.2  (( ff[0] is the occupancy of dorsal )
tsites.clear();
for ( int k = 0; k < motifs.size(); k++ ) {
    for( int i = 0; i < sitesp[k].size(); i++){
        tsites.push_back( sitesp[k][i] );
    }
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
return sitest.size();
}// if snail 0
}// else
int SeqAnnotator::annoty3( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites)
{
    int ncol =e.nCols();
    int i=ti;
    vector< double > reD;
    reD = e.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    vector< double > cons = f.getCol( column );
    tsites.clear();
    vector< double > ff(3);  // try this without initializing
    typedef  vector< Site > sitestt;
    vector< sitestt > sitest( motifs.size() );    // sitestemp
    vector< sitestt > sitesp( motifs.size() );    // sitespermanent
    { int k =2;
    for ( int i = 0; i < seq.size(); i++ ) {
        int l = motifs[ k ].length();
        if ( i + l > seq.size() ) break ;
        double energy;
        Sequence elem( seq, i, l, 1 );
        energy = motifs[ k ].energy( elem );
        if ( energy <= energyThrs[ k ] ) {
            sitest[k].push_back( Site( i, 1, k, energy ) );
        }
        Sequence rcElem( seq, i, l, 0 );
        energy = motifs[ k ].energy( rcElem );
        if ( energy <= energyThrs[ k ] ) {
            sitest[k].push_back( Site( i, 0, k, energy ) );
        }
    }
    double mine =100;
    int gbest;
    int counter=1;
    sitestt  s = sitest[k];
    sitestt  sp;
    sp.clear();
    Site sns;
    for( int i =0; i < seqSites.size() ; i++ ) {
        sns = seqSites[i];
        if( sns.factorIdx == 2 ) {
            sp.push_back( sns );
        }
    }
    sitestoverlap( sitest[2], sp );
    sitestoverlap( s, sp );
    sp.clear();
    while( counter > 0 && s.size() > 0 ){
        counter=0;
        Site tempsite;
        for ( int i = 0; i < s.size(); i++ ) {
            if ( mine >= s[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
            mine = s[i].energy;
            tempsite = s[i];
            gbest = i;
            counter++;
        }
    }
    if( counter != 0) {
        sp.push_back (tempsite);
        sitestoverlap( s, sp );
    }
} // while
s.clear();
for( int i=0; i< seqSites.size(); i++){
    tsites.push_back( seqSites[i] );
}
for( int i=0; i< sp.size(); i++){
    tsites.push_back( sp[i] );
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
double	pre=p.predictExpr( tsites, 5, cons) ;
if(e(ti,column) -e(ti,mi) == 0 ) {sp.clear(); return 1;} //|| e(ti,column)-e(ti,column+2)
if( pre < .5 )
{sp.clear(); return 1; }
else{   sitestoverlap( sitest[2],sp ) ; sp.clear(); } // tsites.clear();}
} // k
int ktempindex=0;
double Etemp=0;
double diffE= 3;
double diffbest=3;
int  kbest=0;
for ( int k = 0; k < sitest[2].size(); k++ ) {
    SiteVec sitesorderedtemp = tsites;  // if this is seqSites, we need to erase tsites before this point, and not kill any of sitest[2]
    ktempindex = k;
    sitesorderedtemp.push_back( sitest[2][k] );
    std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
    Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
    sitesorderedtemp.clear();
    ktempindex=k;
    diffE = abs( Etemp- 0 ) ;    //e(ti,column) );
    if(diffbest > diffE ) {
        diffbest = diffE;
        kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
    }// if
}
sitesp[2].clear();
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] );
int N=1;
double previous = 0;
double current = 0;
if( !sitest[2].empty())  {
    do {
        if( N == 1 ) {
            { int k = 2;
            for( int i = 0; i < sitesp[k].size(); i++){
                tsites.push_back( sitesp[k][i] );
            }
        }
    }// if N
    N++;
    std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
    previous=p.predictExpr( tsites, 5, cons) ;
    if( previous < .5 )
    { break ; }
    if( previous < .2 )
    { break ; }
    if( !sitest[2].empty() ){
        for ( int k = 0; k < sitest[2].size(); k++ ) {
            SiteVec sitesorderedtemp = tsites;
            ktempindex = k;
            sitesorderedtemp.push_back( sitest[2][k] );
            std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
            Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
            sitesorderedtemp.clear();
            ktempindex=k;
            diffE = abs( Etemp-  0 );
            if(diffbest > diffE ) {
                diffbest = diffE;
                kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
            }// if
        }
        tsites.push_back( sitest[2][kbest] ) ;
        std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
        current = p.predictExpr(tsites,5,  cons);
        if( current < .5 ) break;
        sitesp[2].clear();
        sitesp[2].push_back( sitest[2][kbest] );
        sitestoverlap( sitest[2], sitesp[2] );
    } // if sitest[2] not empty
    else break;
} while( N < 3  ); // !sitest[2].empty() );   //abs(previous - current ) > .25 );
} // if sitest2 empty
return tsites.size();
}
int SeqAnnotator::annotyd( const Sequence& seq, SiteVec& tsites, Matrix f, Matrix e, ExprFunc& p, int column, string  names, int ti ,SiteVec& seqSites, vector< SiteVec >& dd)
{
    vector< double > cons = f.getCol( column );
    tsites.clear();
    vector< double > ff(3);  // try this without initializing
    typedef  vector< Site > sitestt;
    vector< sitestt > sitest( motifs.size() );    // sitestemp
    vector< sitestt > sitesp( motifs.size() );    // sitespermanent
    { int k =2;
    for ( int i = 0; i < seq.size(); i++ ) {
        int l = motifs[ k ].length();
        if ( i + l > seq.size() ) break ;
        double energy;
        Sequence elem( seq, i, l, 1 );
        energy = motifs[ k ].energy( elem );
        sitest[k].push_back( Site( i, 1, k, energy ) );
        Sequence rcElem( seq, i, l, 0 );
        energy = motifs[ k ].energy( rcElem );
        sitest[k].push_back( Site( i, 0, k, energy ) );
    }
}
int ktempindex=0;
double Etemp=0;
double diffE= 3;
double diffbest=3;
int  kbest=0;
for ( int k = 0; k < sitest[2].size(); k++ ) {
    SiteVec sitesorderedtemp = seqSites;
    ktempindex = k;
    sitesorderedtemp.push_back( sitest[2][k] );
    std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
    Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
    sitesorderedtemp.clear();
    ktempindex=k;
    diffE = abs( Etemp- 0 ) ;    //e(ti,column) );
    if(diffbest > diffE ) {
        diffbest = diffE;
        kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
    }// if
}
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] );
int N=1;
double previous = 0;
double current = 0;
do {
    tsites.clear();
    tsites = seqSites;
    { int k = 2;
    for( int i = 0; i < sitesp[k].size(); i++){
        tsites.push_back( sitesp[k][i] );
    }
}
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);  // the for loop with sites and this sort may no longer be needed
previous=p.predictExpr( tsites, 5, cons) ;
if( abs(previous - e(ti,column) ) < .2 )
{ break ; }
if( previous < .2 )
{ break ; }
N++;
for ( int k = 0; k < sitest[2].size(); k++ ) {
    SiteVec sitesorderedtemp = tsites;
    ktempindex = k;
    sitesorderedtemp.push_back( sitest[2][k] );
    std::sort(sitesorderedtemp.begin(), sitesorderedtemp.end(), siteSortPredicate);
    Etemp= p.predictExpr( sitesorderedtemp, 5, cons) ;
    sitesorderedtemp.clear();
    ktempindex=k;
    diffE = abs( Etemp- 0 ) ;// e(ti,column) );
    if(diffbest > diffE ) {
        diffbest = diffE;
        kbest = ktempindex;			// this code doesn't need all these variables/ but this way it's easier to read;
    }// if
}
tsites.push_back( sitest[2][kbest] ) ;
std::sort(tsites.begin(), tsites.end(), siteSortPredicate);
current = p.predictExpr(tsites,5,  cons);
if( abs(current - e(ti,column) ) < .2 ) break;
if( current  < .2 ) break;
sitesp[2].push_back( sitest[2][kbest] );
sitestoverlap( sitest[2], sitesp[2] );
} while( N < 3 || abs(previous - current ) > .25 );
Delete1( tsites,dd );  // this is for robustness analysis, checking to see how much the borders shift
return tsites.size();
}
bool SeqAnnotator::DeltaZ( vector< vector< Site > >& sitespa , vector< double >& cons , double Zbefore, ExprFunc& p, double Obefore)
{
    vector< double > f(0);
    vector< Site > tsitess(0);
    tsitess.clear();
    sitespa[1].pop_back();  // remove last twist site
    for ( int k = 0; k < motifs.size(); k++ ) {
        for( int i = 0; i < sitespa[k].size(); i++){
            tsitess.push_back( sitespa[k][i] );
        }
    }
    std::sort(tsitess.begin(), tsitess.end(), siteSortPredicate);
    double Zafter = p.predictZ( tsitess, cons);
    p.predictOcc(tsitess,5,cons,f);
    double Oafter = f[0];
    if ( (Zbefore - Zafter)/Zbefore < .1  && (Obefore - Oafter)/Obefore < .1) { sitespa.clear(); return true; }
    sitespa.clear();
    return false;
}
int SeqAnnotator::Delete1( vector< Site >& t,vector< vector< Site > >& tt )
{
    int size = t.size();
    int count=0;
    for( int i = 0; i < size; i++){
        SiteVec temp1;
        temp1=t;
        vector< Site >::iterator ptr2 = temp1.begin();
        ptr2 += i;
        temp1.erase(ptr2  );
        tt.push_back( temp1 );
    }
    return 1;
}
int SeqAnnotator::maxZ( vector< double >& Z )
{
    int tempindex = 0 ;
    double maxZ=(*min_element(Z.begin(), Z.end() ) );
    for ( int i = 0; i < Z.size(); i++ ) {
        if ( maxZ < Z[i] ){
            maxZ = Z[i];
            tempindex = i;
        }
    }
    return tempindex;
}
double SeqAnnotator::maxZs( vector< double >& Z )
{
    int tempindex = 0 ;
    double maxZ=(*min_element(Z.begin(), Z.end() ) );
    for ( int i = 0; i < Z.size(); i++ ) {
        if ( maxZ < Z[i] ){
            maxZ = Z[i];
            tempindex = i;
        }
    }
    return maxZ;
}
void SeqAnnotator::sitestoverlap( vector< Site >& sitest, vector< Site >& sitesp )
{
    for ( int i = 0; i < sitesp.size(); i++ ) {
        vector< Site >::iterator ptr2 = sitest.begin();
        int j = 0;
        while( ptr2 != sitest.end() ){
            if (siteOverlap( sitesp[i], sitest[j], this->motifs) )  { // overlap returns true when they DO overlap
            sitest.erase(ptr2);
            continue;  // don't increment j, if we erase the site at position j, because the sitest subvector shifts such that a new site fills site j.
        }  // if j = 0, has a 10 bp overlapping site, then we don't want to increment j until the ptr2 loop executes 10 times, for each iteration j remains 0.
        j++;
        ptr2++;
    }
}
}
Site SeqAnnotator::siteMax( vector< Site >& sitess)
{
    Site tempsite;
    double mine=100000;
    for ( int i = 0; i < sitess.size(); i++ ) {
        if ( mine > sitess[i].energy ){  // energy = abs(motif - maxmotif),  (site max looks for lowest energy)
        mine = sitess[i].energy;
        tempsite = sitess[i];
    }
}
return tempsite;
}
int SeqAnnotator::compEnergy( const Sequence& seq, SiteVec& sites ) const
{
    for ( int i = 0; i < sites.size(); i++ ) {
        Sequence elem( seq, sites[i].start, motifs[sites[i].factorIdx].length(), sites[i].strand );
        sites[i].energy = motifs[sites[i].factorIdx].energy( elem );
        sites[i].wtRatio = exp( sites[i].energy-10 );  // 1115 changed - to + in exponential
    }
    return 0;
}
ModelType getModelOption( const string& modelOptionStr )
{
    if ( toupperStr( modelOptionStr ) == "LOGISTIC" ) return LOGISTIC;
    if ( toupperStr( modelOptionStr ) == "DIRECT" ) return DIRECT;
    if ( toupperStr( modelOptionStr ) == "QUENCHING" ) return QUENCHING;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_UNLIMITED" ) return CHRMOD_UNLIMITED;
    if ( toupperStr( modelOptionStr ) == "CHRMOD_LIMITED" ) return CHRMOD_LIMITED;
    if ( toupperStr ( modelOptionStr ) == "BINS" ) return BINS;
    cerr << "modelOptionStr is not a valid model option" << endl;
    exit(1);
}
string getModelOptionStr( ModelType modelOption )
{
    if ( modelOption == LOGISTIC ) return "Logisitic";
    if ( modelOption == DIRECT ) return "Direct";
    if ( modelOption == QUENCHING ) return "Quenching";
    if ( modelOption == CHRMOD_UNLIMITED ) return "ChrMod_Unlimited";
    if ( modelOption == CHRMOD_LIMITED ) return "ChrMod_Limited";
    if ( modelOption == BINS ) return "Bins";
    return "Invalid";
}
string getIntOptionStr( FactorIntType intOption )
{
    if ( intOption == BINARY ) return "Binary";
    if ( intOption == GAUSSIAN ) return "Gaussian";
    if ( intOption == BINSF ) return "Binsf";
    return "Invalid";
}
ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;
    cerr << "objOptionStr is not a valid option of objective function" << endl;
    exit(1);
}
string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";
    return "Invalid";
}
string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";
    return "Invalid";
}
double FactorIntFuncBinsf::compFactorInt( double normalInt, double dist, bool orientation )  const
{
    assert( dist >= 0 );
    return normalInt;
}
double FactorIntFuncBinary::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
    double spacingTerm = ( dist < distThr ? normalInt : 1.0 );
    double orientationTerm = orientation ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
}
double FactorIntFuncGaussian::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
    double GaussianInt = dist < distThr ? normalInt * exp( - ( dist * dist ) / ( 2.0 * sigma * sigma ) ) : 1.0;
    return max( 1.0, GaussianInt );
}
double FactorIntFuncGeometric::compFactorInt( double normalInt, double dist, bool orientation ) const
{
    assert( dist >= 0 );
    double spacingTerm = max( 1.0, dist <= distThr ? normalInt : normalInt * pow( spacingEffect, dist - distThr ) );
    double orientationTerm = orientation ? 1.0 : orientationEffect;
    return spacingTerm * orientationTerm;
}
ExprPar::ExprPar( int _nFactors, vector< vector< vector<double> > > _theV , vector< vector< vector<double> > > _theVr): factorIntMat() , theV(_theV),theVr(_theVr)  // call to Matrix constructor, in the initialization list.  neet to assert that _theV has dimensions of nfactors..
{
    assert( _nFactors > 0 );
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );
    }
    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );
}
ExprPar::ExprPar( int _nFactors ) : factorIntMat()
{
    assert( _nFactors > 0 );
    for ( int i = 0; i < _nFactors; i++ ) {
        maxBindingWts.push_back( ExprPar::default_weight );
    }
    factorIntMat.setDimensions( _nFactors, _nFactors );
    factorIntMat.setAll( ExprPar::default_interaction );
    if (modelOption == BINS ) {theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>(ExprPar::nbins,1)));
    theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>(ExprPar::nbins,1)));}
    for ( int i = 0; i < _nFactors; i++ ) {
        double defaultEffect = modelOption == LOGISTIC ? ExprPar::default_effect_Logistic : ExprPar::default_effect_Thermo;
        txpEffects.push_back( defaultEffect );
    }
}
ExprPar::ExprPar( const vector< double >& _maxBindingWts, const Matrix& _factorIntMat, const vector< double >& _txpEffects, const vector< double >& _repEffects, double _basalTxp ) : maxBindingWts( _maxBindingWts ), factorIntMat( _factorIntMat ), txpEffects( _txpEffects ), repEffects( _repEffects ), basalTxp( _basalTxp )
{
    if ( !factorIntMat.isEmpty() ) assert( factorIntMat.nRows() == maxBindingWts.size() && factorIntMat.isSquare() );
}
ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) : factorIntMat()
{
    int _nFactors = actIndicators.size();
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
    int counter = 0;
    if ( estBindingOption ) {
        for ( int i = 0; i < _nFactors; i++ ) {
            double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) : exp( pars[counter++] );
            maxBindingWts.push_back( weight );
        }
    } else {
        for ( int i = 0; i < _nFactors; i++ ) maxBindingWts.push_back( ExprPar::default_weight );
    }
    if (modelOption == BINS ) {
        theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));
        for ( int i = 0; i < _nFactors; i++ ) {
            for ( int j = 0; j < i; j++ ) {         // 1109 this was changed from <= to <  due to abortion errors, this menas the offdiags are define below
            if ( coopMat( i, j ) ) {
                for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
                double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : exp( pars[counter++] );  //524 restet interaction, was 1, same for factorIntMat
                theV[i][j][k] = interaction;
            }
        }  //if construct
        else {
            for(int k=0; k< ExprPar::nbins; k++){
                theV[i][j][k] = 1;
            }
        }
    }
}
for ( int i = 0; i < _nFactors; i++ ) {
    for ( int j = i ; j < _nFactors; j++ ) {            // this changed from i+1 to i, as indicated in 1109 comment above for offdiags
    for ( int k = 0; k <  ExprPar::nbins; k++ ) {
        theV[i][j][k] =  theV[j][i][k] ;
    }
}
}
theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));
for ( int i = 0; i < _nFactors; i++ ) {
    for ( int j = 0; j < i; j++ ) {
        if ( repIndicators[ i] ) {
            for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
            double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interactionr ), log( max_interactionr ) ) ) : exp( pars[counter++] );  //524 reste]t interaction, was 1, same for factorIntMat
            theVr[i][j][k] = interaction;
        }
    }  //if construct
    else {
        for(int k=0; k< ExprPar::nbins; k++){
            theVr[i][j][k] = 1;
        }
    }
}
}
for ( int i = 0; i < _nFactors; i++ ) {
    for ( int j = i + 1; j < _nFactors; j++ ) {                      // how do we define the diagonal elements, we may need to put an if statement here
    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		// and also start j at i, not i+1, where if(j ==i) then Vr = 1, else vr(ijk) = vr(jik)
    theVr[i][j][k] =  theVr[j][i][k] ;
}
}
}
for ( int i = 0; i < _nFactors; i++ ) {
    double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : exp( pars[counter++] );
    txpEffects.push_back( effect );
}
}
}
ExprPar::ExprPar( const vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators , bool a) : factorIntMat()
{
    int _nFactors = actIndicators.size();
    assert( coopMat.isSquare() && coopMat.nRows() == _nFactors );
    assert( repIndicators.size() == _nFactors );
    int counter = 0;
    if ( estBindingOption ) {
        for ( int i = 0; i < _nFactors; i++ ) {
            double weight = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_weight ), log( max_weight ) ) ) :  pars[counter++] ;
            maxBindingWts.push_back( weight );
        }
    } else {
        for ( int i = 0; i < _nFactors; i++ ) maxBindingWts.push_back( ExprPar::default_weight );
    }
    if (modelOption == BINS ) {
        theV = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));
        for ( int i = 0; i < _nFactors; i++ ) {
            for ( int j = 0; j < i; j++ ) {         // 1109 this was changed from <= to <  due to abortion errors, this menas the offdiags are define below
            if ( coopMat( i, j ) ) {
                for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
                double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interaction ), log( max_interaction ) ) ) : pars[counter++] ;  //524 restet interaction, was 1, same for factorIntMat
                theV[i][j][k] = interaction;
            }
        }  //if construct
        else {
            for(int k=0; k< ExprPar::nbins; k++){
                theV[i][j][k] = 1;
            }
        }
    }
}
for ( int i = 0; i < _nFactors; i++ ) {
    for ( int j = i ; j < _nFactors; j++ ) {            // this changed from i+1 to i, as indicated in 1109 comment above for offdiags
    for ( int k = 0; k <  ExprPar::nbins; k++ ) {
        theV[i][j][k] =  theV[j][i][k] ;
    }
}
}
theVr = vector< vector< vector<double> > >(_nFactors,vector< vector<double> >(_nFactors, vector< double>( ExprPar::nbins,1)));
for ( int i = 0; i < _nFactors; i++ ) {
    for ( int j = 0; j < i; j++ ) {
        if ( repIndicators[ i] ) {
            for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
            double interaction = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_interactionr ), log( max_interactionr ) ) ) : pars[counter++] ;  //524 reste]t interaction, was 1, same for factorIntMat
            theVr[i][j][k] = interaction;
        }
    }  //if construct
    else {
        for(int k=0; k< ExprPar::nbins; k++){
            theVr[i][j][k] = 1;
        }
    }
}
}
for ( int i = 0; i < _nFactors; i++ ) {
    for ( int j = i + 1; j < _nFactors; j++ ) {                      // how do we define the diagonal elements, we may need to put an if statement here
    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		// and also start j at i, not i+1, where if(j ==i) then Vr = 1, else vr(ijk) = vr(jik)
    theVr[i][j][k] =  theVr[j][i][k] ;
}
}
}
for ( int i = 0; i < _nFactors; i++ ) {
    double effect = searchOption == CONSTRAINED ? exp( inverse_infty_transform( pars[counter++], log( min_effect_Thermo ), log( max_effect_Thermo ) ) ) : pars[counter++] ;
    txpEffects.push_back( effect );
}
}
}
void ExprPar::getFreePars( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            pars.push_back( weight );
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
                for(int k=0; k< ExprPar::nbins ; k++){
                    double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) : log( theV[i][j][k] );
                    pars.push_back( interaction);
                }
            }  //if construct
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j < i; j++ ) {   // 1026 this was change from <= to < since too many pars were found
        if ( repIndicators[i] ) {
            for(int k=0; k< ExprPar::nbins ; k++){
                double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : log( theVr[i][j][k] );
                pars.push_back( interaction);
            }
        }  //if construct
    }
}
for ( int i = 0; i < nFactors(); i++ ) {
    double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : log( txpEffects[i] );
    pars.push_back( effect );
}
}
void ExprPar::getFreePars3( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators ) const
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : maxBindingWts[i] ;
            pars.push_back( weight );
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
                for(int k=0; k< ExprPar::nbins ; k++){
                    double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) :  theV[i][j][k] ;
                    pars.push_back( interaction);
                }
            }  //if construct
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j < i; j++ ) {   // 1026 this was change from <= to < since too many pars were found
        if ( repIndicators[i] ) {
            for(int k=0; k< ExprPar::nbins ; k++){
                double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : theVr[i][j][k] ;
                pars.push_back( interaction);
            }
        }  //if construct
    }
}
for ( int i = 0; i < nFactors(); i++ ) {
    double effect = searchOption == CONSTRAINED ? infty_transform( log( txpEffects[i] ), log( min_effect_Thermo ), log( max_effect_Thermo ) ) : txpEffects[i] ;
    pars.push_back( effect );
}
}
void ExprPar::getFreePars2( vector< double >& pars, const IntMatrix& coopMat, const vector< bool >& actIndicators, const vector< bool >& repIndicators )
{
    assert( coopMat.isSquare() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() && repIndicators.size() == nFactors() );
    pars.clear();
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            double weight = searchOption == CONSTRAINED ? infty_transform( log( maxBindingWts[ i ] ), log( min_weight ), log( max_weight ) ) : log( maxBindingWts[i] );
            pars.push_back( weight );
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if (coopMat(i,j) ) {
                for(int k=0; k< ExprPar::nbins ; k++){
                    double interaction = searchOption == CONSTRAINED ? infty_transform( log( theV[i][j][k] ), log( min_interaction ), log( max_interaction ) ) : log( theV[i][j][k] );
                    pars.push_back( interaction);
                }
            }  //if construct
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j < i; j++ ) {   // 1026 this was change from <= to < since too many pars were found
        if ( repIndicators[i] ) {
            for(int k=0; k< ExprPar::nbins ; k++){
                double interaction = searchOption == CONSTRAINED ? infty_transform( log( theVr[i][j][k] ), log( min_interactionr ), log( max_interactionr ) ) : log( theVr[i][j][k] );
                pars.push_back( interaction);
            }
        }  //if construct
    }
}
for ( int i = 0; i < nFactors(); i++ ) {
    double effect = searchOption == CONSTRAINED ? infty_transform(  txpEffects[i] , min_effect_Thermo ,  max_effect_Thermo  ) :  txpEffects[i] ;
    pars.push_back( effect );
}
}
void ExprPar::print( ostream& os, const vector< string >& motifNames, const IntMatrix& coopMat ) const //, const vector< vector< vector<double > > > _theV ) const
{
    for ( int i = 0; i < nFactors(); i++ ) {
        os << "motif"<<motifNames[i] << "\t" << "mBwt" << maxBindingWts[i] << "\t" <<"txpE"<< txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) os << "\t" << repEffects[i];
        os << endl;
    }
    if(modelOption== BINS){
        for(int i =0; i < theV.size(); i++){
            for(int j =0; j < theV[i].size(); j++){
                for(int k =0; k < theV[i][j].size(); k++)  //  i don't think this should be k less then vij.size...
                {  if ( coopMat(i,j)){
                    os << "theV" << "\t" << theV[i][j][k] ;
                }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
            }
        }
    }
}
}
int ExprPar::load( const string& file, IntMatrix& coopMat, const vector< bool >& repIndicators )
{
    ifstream fin( file.c_str() );
    if ( !fin ){ cerr << "Cannot open parameter file " << file << endl;	exit( 1 ); }
    vector< string > motifNames( nFactors() );
    for ( int i = 0; i < nFactors(); i++ ) {
        fin >> motifNames[i] >> maxBindingWts[i] >> txpEffects[i];
        if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) fin >> repEffects[i];
    }
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
    theV = vector< vector< vector<double> > >(nFactors(),vector< vector<double> >(nFactors(), vector< double>( ExprPar::nbins,1)));
    theVr = vector< vector< vector<double> > >(nFactors(),vector< vector<double> >(nFactors(), vector< double>( ExprPar::nbins,1)));
    string value, sym;
    fin >> sym;
    fin >> sym;
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) {
                for(int k=0; k< ExprPar::nbins; k++){
                    fin >> value;
                    double interaction = atof( value.c_str() ) ; // for more complicated data use factorIdxMap[motifNames[i]] for i,j indices inside for loops
                    theV[i][j][k] = interaction;
                }
            }  //if construct
            else {
                for(int k=0; k< ExprPar::nbins; k++){
                    theV[i][j][k] = 1;
                }
            }
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = i + 1; j < nFactors(); j++ ) {
            for ( int k = 0; k <  ExprPar::nbins; k++ ) {
                theV[i][j][k] =  theV[j][i][k] ;
            }
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( repIndicators[ i] ) {
                fin >> sym;
                fin >> sym;
                for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
                fin >> value;
                double interaction = atof( value.c_str() ) ; // for more complicated data use factorIdxMap[motifNames[i]] for i,j indices inside for loops
                theVr[i][j][k] = interaction;
            }
        }  //if construct
        else {
            for(int k=0; k< ExprPar::nbins; k++){
                theVr[i][j][k] = 1;
            }
        }
    }
}
for ( int i = 0; i < nFactors(); i++ ) {
    for ( int j = i + 1; j < nFactors(); j++ ) {                      // how do we define the diagonal elements, we may need to put an if statement here
    for ( int k = 0; k <  ExprPar::nbins; k++ ) {		// and also start j at i, not i+1, where if(j ==i) then Vr = 1, else vr(ijk) = vr(jik)
    theVr[i][j][k] =  theVr[j][i][k] ;
}
}
}
return 0;
}
void ExprPar::adjust()
{
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) ) maxBindingWts[i] *= 2.0;
        if ( maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) ) maxBindingWts[i] /= 2.0;
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            for(int k =0; k < ExprPar::nbins; k++){
                if ( theV[i][j][k] < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) ) {
                    theV[i][j][k] *= 10.0;
                    theV[i][j][k] = theV[j][i][k] ;
                }
                if ( theV[i][j][k] > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) {
                    theV[i][j][k] /= 2.0;
                    theV[i][j][k] = theV[j][i][k] ;
                }
            }
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            for(int k =0; k < ExprPar::nbins; k++){
                if ( theVr[i][j][k] < ExprPar::min_interactionr * ( 1.0 + ExprPar::delta ) ) {
                    theVr[i][j][k] *= 10.0;
                    theVr[i][j][k] = theVr[j][i][k] ;
                }
                if ( theVr[i][j][k] > ExprPar::max_interactionr * ( 1.0 - ExprPar::delta ) ) {
                    theVr[i][j][k] /= 2.0;
                    theVr[i][j][k] = theVr[j][i][k] ;
                }
            }
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) ) txpEffects[i] *= 2.0;
        if ( txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) ) txpEffects[i] /= 2.0;
    }
}
ModelType ExprPar::modelOption = CHRMOD_UNLIMITED;  //   5/19
SearchType ExprPar::searchOption = UNCONSTRAINED;   //   5/19
int ExprPar::estBindingOption = 1;                  //   1. estimate binding parameters; 0. not estimate binding parameters
vector< int > ExprFunc::Bc(5,1);
vector< Sequence > ExprPredictor::seqsy;
vector< string > ExprPredictor::seqNmes;
Matrix ExprPredictor::factorExprData2;
Matrix ExprPredictor::exprData2;
double ExprPar::default_weight = 10;
double ExprPar::default_interaction = 1;
double ExprPar::default_effect_Logistic = 5;
double ExprPar::default_effect_Thermo = 5;
double ExprPar::default_repression = 1.0E-2;
double ExprPar::default_basal_Logistic = -9.0;
double ExprPar::default_basal_Thermo = 0.0001;
double ExprPar::min_weight = .99999;
double ExprPar::max_weight = 100;
double ExprPar::min_interaction = .99999;
double ExprPar::max_interaction = 100;
double ExprPar::min_interactionr = .001;
double ExprPar::max_interactionr = 1.00001;
double ExprPar::min_effect_Logistic = 5;
double ExprPar::max_effect_Logistic = 5;
double ExprPar::min_effect_Thermo = 1; //.01;
double ExprPar::max_effect_Thermo = 100;
double ExprPar::min_repression = 1.0E-3;
double ExprPar::max_repression = 500;
double ExprPar::min_basal_Logistic = -9.0;
double ExprPar::max_basal_Logistic = -9.0;
double ExprPar::min_basal_Thermo =.001;
double ExprPar::max_basal_Thermo = 0.001;
double ExprPar::delta = 0.0000001;  // was .1 1010
int ExprPar::nbins = 5;  //binwidth;
ExprFunc::ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par,  const vector < int >& _B, const vector < int >& _Br  ) : motifs( _motifs ), intFunc( _intFunc ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), par( _par ),   B( _B ), Br( _Br )
{
    int nFactors = par.nFactors();
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );
    Bc = vector< int >(5,1);
}
void ExprFunc::predictOcc( const SiteVec& _sites, int length, const vector< double >& factorConcs , vector< double >& fOcc)  // k is the current sequence
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j;
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break;
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[ sites[i].factorIdx ] * sites[i].wtRatio );
    }
    if ( modelOption == BINS ) {
        vector< double > Y( n + 1 );
        Y[0] = 0;
        vector< double > Yt( n + 1 );
        Yt[0] = 0;
        double Z_bind = compPartFunc();
        for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
        for ( int i = 1; i <= n; i++ ) {
            int factorIdx;
            factorIdx = k;
            int match = sites[ i ].factorIdx == factorIdx ? 1 : 0;
            double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
            for ( int j = boundaries[i] + 1; j < i; j++ ) {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
            }
            Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
            Yt[i] = Y[i] + Yt[i - 1];
        }
        double Y_total = Yt[n];
        factorOcc.push_back( Y_total / Z_bind );
    }
    double totalEffect = 0;
    fOcc = factorOcc;
    for ( int i = 0; i < motifs.size(); i++ ) {
        double effect = factorOcc[i];  // par.txpEffects[i] * factorOcc[i]; // 103  // we need to make sure the order of par.txpEffects and factorOcc are the same for the factors, i.e i=i?
        totalEffect += effect;
    }
}
}
double ExprFunc::predictExpr2( const SiteVec& _sites, int length, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j;
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break;
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );
    }
    vector< double > Y( n + 1 );
    Y[0] = 0;
    vector< double > Yt( n + 1 );
    Yt[0] = 0;
    double Z_bind = compPartFunc();
    for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
    for ( int i = 1; i <= n; i++ ) {
        int factorIdx;
        factorIdx = k;
        int match = sites[ i ].factorIdx == factorIdx ? 1 : 0;
        double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
        }
        Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
        Yt[i] = Y[i] + Yt[i - 1];
    }
    double Y_total = Yt[n];
    factorOcc.push_back( Y_total / Z_bind );
}
double totalEffect = 0;
totalEffect = factorOcc[0]  ;// -5;
if(totalEffect >= .5){  return  totalEffect; }
else{ return 0;}
}
double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j;
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break;
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );
    }
    vector< double > Y( n + 1 );
    Y[0] = 0;
    vector< double > Yt( n + 1 );
    Yt[0] = 0;
    double Z_bind = compPartFunc();
    for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
    for ( int i = 1; i <= n; i++ ) {
        int factorIdx;
        factorIdx = k;
        int match = sites[ i ].factorIdx == factorIdx ? 1 : 0;
        double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
        }
        Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
        Yt[i] = Y[i] + Yt[i - 1];
    }
    double Y_total = Yt[n];
    factorOcc.push_back( Y_total / Z_bind );
}
double totalEffect = 0;
totalEffect = factorOcc[0]  ;// -5;
return  totalEffect;  //  1/(1+exp(-totalEffect) );                    //logistic(  totalEffect-1 );
}
double ExprFunc::predictExpr3( const SiteVec& _sites, int length, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j;
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break;
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );
    }
    vector< double > Y( n + 1 );
    Y[0] = 0;
    vector< double > Yt( n + 1 );
    Yt[0] = 0;
    double Z_bind = compPartFunc();
    for ( int k = 0; k < motifs.size(); k++ ) {  // should k start with 0 or 1?  101jc  (we know what k is by the order it was placed in the file  )
    for ( int i = 1; i <= n; i++ ) {
        int factorIdx;
        factorIdx = k;
        int match = sites[ i ].factorIdx == factorIdx ? 1 : 0;
        double sum = Yt[boundaries[i]] + match * Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * ( Y[ j ] +  match * Z[ j ] );
        }
        Y[ i ] = bindingWts[ i ] * sum;    // here is where the concentrations are being accounted for factor i.  what about for factor j, those too, cr bindingWts.
        Yt[i] = Y[i] + Yt[i - 1];
    }
    double Y_total = Yt[n];
    factorOcc.push_back( Y_total / Z_bind );
}
double totalEffect = 0;
totalEffect =par.txpEffects[0] * factorOcc[0] + par.txpEffects[1]*factorOcc[1]  -  par.txpEffects[2] * factorOcc[2] -5;
return  1/(1+exp(-totalEffect) );                    //logistic(  totalEffect-1 );
}
double ExprFunc::predictZ( const SiteVec& _sites, const vector< double >& factorConcs )
{
    bindingWts.clear();
    boundaries.clear();
    factorOcc.clear();
    int n = _sites.size();
    sites = _sites;
    sites.insert( sites.begin(), Site() );  // start with a pseudo-site at position 0
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ ) {
        int j;
        for ( j = i - 1; j >= 1; j-- ) {
            if ( ( sites[i].start - sites[j].start ) > range ) break;
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ ) {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].wtRatio );
    }
    double Z_bind = compPartFunc();
    return Z_bind;
}
double ExprFunc::compPartFunc()
{
    int n = sites.size() - 1;
    Z = vector< double >( n + 1 );
    Z[0] = 1.0;
    Zt = vector< double >( n + 1 );
    Zt[0] = 1.0;
    double sum=0;  //1109
    for ( int i = 1; i <= n; i++ ) {
        sum = Zt[boundaries[i]];  //  declaration for sum was here inside for loop
        for ( int j = boundaries[i] + 1; j < i; j++ ) {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }
        Z[i] = bindingWts[ i ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }
    double Z_bind = Zt[n];
    return Z_bind;
}
ModelType ExprFunc::modelOption = QUENCHING;  //5/19
double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
    double dist = abs( a.start - b.start );
    double maxInt=1;
    if ( ( repIndicators[a.factorIdx ] || repIndicators[ b.factorIdx ] ) && ( repIndicators[a.factorIdx ] && repIndicators[ b.factorIdx ] )   )
    maxInt = 1;
    else if ( repIndicators[ a.factorIdx ] || repIndicators[ b.factorIdx ]) {
        if (      dist <= Br[0] )                 { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 0 ]; } //both sides of matrix diagonal need  to have the same parameters
        if (      dist > Br[0] && dist <= Br[1] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 1 ]; }
        if (      dist > Br[1] && dist <= Br[2] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 2 ]; }
        if (      dist > Br[2] && dist <= Br[3] ) { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 3 ]; }
        if (      dist > Br[3] )                 { maxInt = par.theVr[a.factorIdx][b.factorIdx][ 4 ]; }
    }
    else {
        if (      dist <= B[0] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 0 ]; ExprFunc::Bc[0]++; }
        if (      dist > B[0] && dist <= B[1] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 1 ]; ExprFunc::Bc[1]++;}
        if (      dist > B[1] && dist <= B[2] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 2 ]; ExprFunc::Bc[2]++;}
        if (      dist > B[2] && dist <= B[3] ) { maxInt = par.theV[a.factorIdx][b.factorIdx][ 3 ]; ExprFunc::Bc[3]++;}
        if (      dist > B[3] )                 { maxInt = par.theV[a.factorIdx][b.factorIdx][ 4 ]; ExprFunc::Bc[4]++;}
    }
    bool orientation = ( a.strand == b.strand );
    return  intFunc->compFactorInt( maxInt, dist, orientation );
}
bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
    double dist = abs( a.start - b.start );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}
ExprPredictor::ExprPredictor(const vector< vector< SiteVec > >& _seqSitesb, const vector< vector< double > >& _bindingData,  vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const FactorIntFunc* _intFunc, const IntMatrix& _coopMat, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr,const vector< int >& _binBord ,const vector< int >& _binBordr, const vector < bool >& _indicator_bool, SeqAnnotator _anny, const string& _file,  vector< SiteVec >& _seqSitesbot,  vector< SiteVec >& _seqSitesm1,  vector< SiteVec >& _seqSitesm2 , vector< SiteVec >& _seqSitesf2 ,vector< SiteVec >& _seqSitesbotf2, vector< SiteVec >& _seqSitesm1f2, vector< SiteVec >& _seqSitesm2f2, vector< SiteVec >& _seqSitesf3 ,  vector< SiteVec >& _seqSitesbotf3,vector< SiteVec >& _seqSitesm1f3 , vector< SiteVec >& _seqSitesm2f3): seqSitesb( _seqSitesb ), bindingData( _bindingData), seqSites( _seqSites ), seqLengths( _seqLengths ), exprData( _exprData ), motifs( _motifs ), factorExprData( _factorExprData ), intFunc( _intFunc ), coopMat( _coopMat ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr ), binBord(_binBord), binBordr(_binBordr),  indicator_bool ( _indicator_bool ), anny( _anny) , file( _file ) , seqSitesbot( _seqSitesbot ), seqSitesm1( _seqSitesm1 ), seqSitesm2( _seqSitesm2 ),seqSitesf2( _seqSitesf2 ) ,seqSitesbotf2(_seqSitesbotf2), seqSitesm1f2( _seqSitesm1f2) ,seqSitesm2f2( _seqSitesm2f2), seqSitesf3( _seqSitesf3) , seqSitesbotf3( _seqSitesbotf3),seqSitesm1f3( _seqSitesm1f3 ), seqSitesm2f3( _seqSitesm2f3), seqSitesm1d1(),d(),spaceweights(0)
{
    assert( exprData.nRows() == nSeqs() );
    assert( factorExprData.nRows() == nFactors() && factorExprData.nCols() == nConds() );
    assert( coopMat.isSquare() && coopMat.isSymmetric() && coopMat.nRows() == nFactors() );
    assert( actIndicators.size() == nFactors() );
    assert( maxContact > 0 );
    assert( repIndicators.size() == nFactors() );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors() );
    assert( repressionDistThr >= 0 );
    ExprPar::modelOption = modelOption;
    ExprFunc::modelOption = modelOption;
    ExprPar::nbins = _binBord.size() +1;  //  add one since there exist a bin after the last border and before the first border
    if ( modelOption != LOGISTIC && modelOption != DIRECT && modelOption !=BINS ) {
        ExprPar::min_effect_Thermo = 0.99;
        ExprPar::min_interaction = 0.99;
    }
    ExprPar::estBindingOption = estBindingOption;
}
double ExprPredictor::objFunc( const ExprPar& par ) const
{
    if ( objOption == SSE ) return compRMSE( par );
    if ( objOption == CORR ) return -compAvgCorr( par );
    if ( objOption == CROSS_CORR ) return -compAvgCrossCorr( par );
}
double ExprPredictor::objFunc2(  ExprPar& par )
{
    return -compAvgCorr2( par );
}
double ExprPredictor::objFuncborder(  ExprPar& par )
{
    return compAvgCorrborder8( par );
}
double ExprPredictor::objFuncborder2(  ExprPar& par )
{
    return compAvgCorrborder2( par );  // changed to + sign 11/6/2011
}
int ExprPredictor::train( const ExprPar& par_init )
{
    par_model = par_init;
    if ( nAlternations == 0 ) {obj_model = objFuncborder2( par_model );return 0;}
    ExprPar par_result;
    double obj_result;
    for ( int i = 0; i < nAlternations; i++ ) {
        gradient_minimize( par_result, obj_result );
        par_model = par_result;
    }
    par_model = par_result;
    obj_model = obj_result;
    printPar( par_model );
    return 0;
}
int ExprPredictor::train4( const ExprPar& par_init )
{
    par_model = par_init;
    cout << " inside train4 " << endl;
    printPar( par_model );
    ExprPar par_result;
    double obj_result;
    gradient_minimize2( par_result, obj_result );
    par_model = par_result;
    obj_model = obj_result;
    printPar( par_model );
    return 0;
}
int ExprPredictor::train( const ExprPar& par_init, const gsl_rng* rng )
{
    train( par_init );
    ExprPar par_best = par_model;
    double obj_best = obj_model;
    for ( int i = 0; i < nRandStarts; i++ ) {
        ExprPar par_curr = par_init;
        randSamplePar( rng, par_curr ); /////////////this has a bug 525
        train( par_curr );
        printPar( par_model );
        if ( obj_model < obj_best ) {
            par_best = par_model;
            obj_best = obj_model;
        }
    }
    if ( nRandStarts ) train( par_best );
    par_model = par_best;
    obj_model = obj_best;
    cout<< " about to print best par " << endl;
    printPar( par_model );
    cout <<" ending train rng " << endl;
    return 0;
}
int ExprPredictor::train3( const ExprPar& par_init, const gsl_rng* rng )
{
    par_model = par_init;
    if ( ExprPar::searchOption == CONSTRAINED ) par_model.adjust();
    obj_model = objFunc( par_model );
    if ( nAlternations == 0 ) return 0;
    return 0;
}
void ExprPredictor::printFile( ostream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    //    cout.precision( 3 );
    //     cout.width( 8 );
    // print binding weights
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        if ( modelOption == LOGISTIC || modelOption == DIRECT ) os << par.txpEffects[i] << "\t";
        else {
            if ( actIndicators[i] ) os << par.txpEffects[i] << "\t";
        }
    }
    // print the basal transcription
    os << par.basalTxp << "\t";
    // print the interaction matrix
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) os << par.factorIntMat( i, j ) << "\t";
        }
    }
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( repIndicators[i] ) os << par.repEffects[i] << "\t";
        }
    }
    os << jRMSE.getObj() << endl;
}
void ExprPredictor::printFilePar_KfoldCV( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 );
    //os.width( 18 );
    // print binding weights
    os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    os << endl<<" txpEffects : "<< endl;
    // print the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ ) {
        os << par.txpEffects[i] << "\t";
    }
    os << endl<<  " coopertivity : " << endl;
    for(int i =0; i < nFactors(); i++){
        for(int j =0; j <= i; j++){
            if ( coopMat(i,j)){
                for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
                os <<  "\t" << par.theV[i][j][k] ;
            }
        }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
        //count = count + 1 ;//ouut << endl;
    }
}
os<< endl << " quenching : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <i; j++){
        if(repIndicators[i])  { cout << endl;
        for(int k =0; k < ExprPar::nbins ; k++){
            os <<  "\t" << par.theVr[i][j][k] ;
        }
    }
}
}
os << endl << "objective function : " << endl;
os << jRMSE.getObj() << endl;
ofstream fo( "oData/format.tex",ios::app );
int j=0;
cout << seqSites.size() << endl;
vector< Site > tsites(0);
fo << endl<<  " cooperativity : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <= i; j++){
        if ( coopMat(i,j)){
            for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
            fo <<  "\t" << par.theV[i][j][k] ;
        }
    }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
    //count = count + 1 ;//ouut << endl;
}
}
fo << endl;
for(int m = 0; m < seqSites.size(); m++ ) {
    tsites = seqSites[m];
    ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
    foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
    for( int i = 0; i < tsites.size() ; i++ ) {
        if (tsites.size() == 1 ) { break; }
        for (;;) {
            if ( j < tsites[i].start ) { foo << "b"; j++;  }
            else {//cout << " j = " << j << endl;
            if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
            if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
            if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
        }// else
    }
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo.close();
/////////////////
j=0;
for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
    if( i == 0 ) { fo << "&&" ; continue; }
    fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
    if(i% 90 == 0 ) {fo << "\\\\&&"; }
}
fo << "\\\\&&\\\\";
fo << ExprPredictor::seqNmes[m] << "&" << "cell"<< cell[m] << "&" ;
for( int i = 0; i < tsites.size() ; i++ ) {
    if (tsites.size() == 1 ) { break; }
    for (;;) {
        if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\&&\\\\&&\\\\" << endl;
}// for m
fo.close();
}
int ExprPredictor::load( const string& fila )
{
    ofstream fout( fila.c_str() );
    vector< string > motifNames;
    string dorsal, twist, snail;
    dorsal="dl";
    twist="tw";
    snail="sn";
    motifNames.push_back( dorsal );
    motifNames.push_back( twist );
    motifNames.push_back( snail );
    fout << setprecision(2) ;
    for ( int i = 0; i < nFactors(); i++ ) {
        fout << motifNames[i] << '\t'<< par_model.maxBindingWts[i]<< '\t' << par_model.txpEffects[i] << '\n';
    }
    map< string, int > factorIdxMap;
    for ( int i = 0; i < nFactors(); i++ ) {
        factorIdxMap[motifNames[i]] = i;
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            if ( coopMat( i, j ) ) {
                fout << motifNames[j]<< '\t' << motifNames[i];
                for(int k=0; k< ExprPar::nbins; k++){
                    fout << '\t' << par_model.theV[i][j][k];
                }
                fout << '\n' ;
            }  //if construct
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j < i; j++ ) {
            if ( repIndicators[ i] ) {
                fout << motifNames[j]<< '\t' << motifNames[i];
                for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
                fout<< '\t'<< par_model.theVr[i][j][k]  ;
            }
            fout << '\n' ;
        }  //if construct
    }
}
fout.close();
return 0;
}
void ExprPredictor::printFiled( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 );
    os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    os << endl<<" txpEffects : "<< endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        os << par.txpEffects[i] << "\t";
    }
    os << endl<<  " coopertivity : " << endl;
    for(int i =0; i < nFactors(); i++){
        for(int j =0; j <= i; j++){
            if ( coopMat(i,j)){
                for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
                os <<  "\t" << par.theV[i][j][k] ;
            }
        }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
    }
}
os<< endl << " quenching : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <i; j++){
        if(repIndicators[i])  { cout << endl;         /// why is this cout here?
        for(int k =0; k < ExprPar::nbins ; k++){
            os <<  "\t" << par.theVr[i][j][k] ;
        }
    }
}
}
os << endl << "objective function : " << endl;
os << jRMSE.getObj() << endl;
ofstream fo( "oData/format.tex",ios::app );
int j=0;
cout << seqSites.size() << endl;
vector< Site > tsites;
vector< Site > tsitesbot;
for(int m = 0; m < seqSites.size(); m++ ) {
    tsites = seqSites[m];
    tsitesbot = seqSitesbot[m]  ;
    ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
    foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
    for( int i = 0; i < tsites.size() ; i++ ) {
        for (;;) {
            if ( j < tsites[i].start ) { foo << "b"; j++;  }
            else { cout << " j1 = " << j << endl;
            if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
            if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
            if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
        }// else
    }
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {
        if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
        else { cout << " j = " << j << endl;
        if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo.close();
j=0;
for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
    if( i == 0 ) { fo << "&&" ; continue; }
    fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
    if(i% 90 == 0 ) {fo << "\\\\&&"; }
}
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
vector< double > concs = factorExprData.getCol( cell[2*m] );
double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< predicted  << "&" << "pcell "<< cell[2*m] << "&" ;
for( int i = 0; i < tsites.size() ; i++ ) {
    for (;;) {
        if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {
        if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 0 << "&" ;
for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 10 << "&" ;
for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
for( int e = 0; e < d[m].size() ; e++ ){  //d[m].size = seqSitesm1[m].size  // since seqSitesm1 is parsed into  all fragments with one missing site
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< 0 << "&" ;
for( int i = 0; i < d[m][e].size() ; i++ ) {
    for (;;) {
        if ( j < d[m][e][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( d[m][e][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( d[m][e][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( d[m][e][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
} // for e
fo << "\\\\&&\\\\&&\\\\" << endl;
}// for m
fo.close();
}
void ExprPredictor::printFile3b( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE )
{
    os.setf( ios::fixed );
    os.precision( 3 );
    os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    os << endl<<" txpEffects : "<< endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        os << par.txpEffects[i] << "\t";
    }
    os << endl<<  " coopertivity : " << endl;
    for(int i =0; i < nFactors(); i++){
        for(int j =0; j <= i; j++){
            if ( coopMat(i,j)){
                for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
                os <<  "\t" << par.theV[i][j][k] ;
            }
        }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
    }
}
os<< endl << " quenching : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <i; j++){
        if(repIndicators[i])  {
            for(int k =0; k < ExprPar::nbins ; k++){
                os <<  "\t" << par.theVr[i][j][k] ;
            }
        }
    }
}
os << endl << "objective function : " << endl;
os << jRMSE.getObj() << endl;
ofstream fo( "oData/format.tex",ios::app );
int j=0;
vector< Site > tsites;
vector< Site > tsitesbot;
double bottomofdborder = .5;
double topofdborder = .8;
int nrow = factorExprData.nRows();
int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
ofstream hes( "hes.txt" );
vector< double > pars;
par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators );
int pars_size = pars.size();
fix_pars.clear();
free_pars.clear();
for( int index = 0; index < pars_size; index++ ){
    if( indicator_bool[ index ] ){
        free_pars.push_back( pars[ index ]);
    }
    else{
        fix_pars.push_back( pars[ index ] );
    }
}
pars.clear();
pars = free_pars;
double step = .01;
bool a=1;
ofstream occm("occmat.txt");
cout << AllData.size() << endl;
cout << "AllData[0].size()"<< AllData[0].size() << endl;
cout << "AllData[1].size()"<< AllData[1].size() << endl;
cout << "AllBorders.size() " << AllBorders.size() << endl;
vector< double > initilizer( pars.size() );
Jacobian = initilizer;
vector< vector< double > > init2( pars.size(), initilizer );
Hessian = init2;
gsl_matrix *JM;  // number of [sequences * (number of bords)] by the number of parameters
int count = 0; // used to indicate data point within bords*seq
JM = gsl_matrix_alloc( AllData[0].size()* AllData.size(),pars.size() );  // number of sequences(*number of bords) by the number of parameters
for(int bords = 0; bords< AllData.size(); bords++) {  //bords represent vectors like siteVec.. each bord has all sequences sequence ( for given bord, loop over all seqs)
for(int i=0; i < AllData[bords].size(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
for (int jj = 0; jj < pars.size() ; jj++ ) {
    vector< double > parsjuh =pars;
    vector< double > parsjdh= pars;
    parsjuh[jj] = pars[jj] + step;
    parsjdh[jj] = pars[jj] - step;
    vector< double > concs = factorExprData.getCol( AllBorders[i][bords] );
    vector < double > all_pars;
    all_pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    for( int index = 0; index < indicator_bool.size(); index ++ ){
        if(  indicator_bool[ index ]  ){
            all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
        }
        else{
            all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
        }
    }
    ExprPar parjuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators , a);
    all_pars.clear();
    free_par_counter = 0;
    fix_par_counter = 0;
    for( int index = 0; index < indicator_bool.size(); index ++ ){
        if(  indicator_bool[ index ]  ){
            all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
        }
        else{
            all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
        }
    }
    ExprPar parjdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
    vector< double > fOcc(0);
    hes.setf( ios::fixed);
    hes.width(3);
    count = bords* nSeqs();
    ExprFunc* func = this->createExprFunc( parjuh );
    double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
    func = this->createExprFunc( parjdh );
    double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
    func = this->createExprFunc( par_model );
    if( jj==0){
        func->predictOcc( AllData[bords][i], 5, concs,  fOcc );
        ;
        occm << ExprPredictor::seqNmes[i]<< '\t'<<AllBorders[i][bords] <<'\t' << fOcc << endl;
    } //if i ==0
    double p = func->predictExpr(AllData[bords][i], 5, concs ); // AllData[bords][i] = seqSites[i]
    Jacobian[jj] += ( pjuh- pjdh )/(2*step);
    gsl_matrix_set( JM,count+ i, jj, ( pjuh- pjdh )/(2*step) );  // i is seqs, jj is parameters
    parsjuh.clear();
    parsjdh.clear();
}  // for jj
for (int j = 0; j < pars.size() ; j++ ) {
    for (int k = 0; k < pars.size() ; k++ ) {
        vector< double > parsjuh =pars;
        vector< double > parsjdh= pars;
        vector< double > parskuh =pars;  // this and the next aren't being used
        vector< double > parskdh= pars;
        vector< double > parsjkdh= pars;
        vector< double > parsjkuh =pars;
        vector< double > parsjukdh= pars;
        vector< double > parsjdkuh =pars;
        parsjuh[j] = pars[j] + step;    // f_j    first partial
        parsjdh[j] = pars[j] - step;    // f_j
        parskuh[k] = pars[k] + step;    // f_k
        parskdh[k] = pars[k] - step;     // f_k
        parsjkuh[j] = pars[j] + step;    //  f_kj   seconder order partial, first step in construction
        parsjkuh[k] = pars[k] + step;     // f_kj  second step in construction: f(x+h,y+h)
        parsjkdh[j] = pars[j] - step;    // f_kj  second order,..  f(k-h,j-h)
        parsjkdh[k] = pars[k] - step;
        parsjukdh[j] = pars[j] + step;     //  f_kj  : f(k-h, j+h)
        parsjukdh[k] = pars[k] - step;
        parsjdkuh[j] = pars[j] - step;     //  f_kj  : f(k+h, j-d )
        parsjdkuh[k] = pars[k] + step;
        vector< double > concs = factorExprData.getCol( AllBorders[i][bords] );
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parjuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parjdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parskuh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parskdh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parkdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        if( j == k ) {
            ExprFunc* func = this->createExprFunc( parjuh );
            double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjdh );
            double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( par_model );
            double p = func->predictExpr(AllData[bords][i], 5, concs );
            Hessian[j][k] += ( pjuh -2*p + pjdh )/(step*step) ;
            parsjuh.clear();
            parsjdh.clear();
            parskuh.clear();
            parskdh.clear();
            parsjkuh.clear();
            parsjkdh.clear();
            parsjdkuh.clear();
            parsjukdh.clear();
        }
        if ( j !=k ) {
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjkuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjkdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjkdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjdkuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjdkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjukdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjukdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            ExprFunc* func = this->createExprFunc( parjuh );
            double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjdh );
            double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( par_model );
            double p = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parkuh );
            double pkuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parkdh );
            double pkdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjkuh );
            double pjkuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjkdh );
            double pjkdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjukdh );
            double pjukdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjdkuh );
            double pjdkuh = func->predictExpr(AllData[bords][i], 5, concs );
            Hessian[j][k] +=( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step);
            parsjuh.clear();
            parsjdh.clear();
            parskuh.clear();
            parskdh.clear();
            parsjkuh.clear();
            parsjkdh.clear();
        }  // if j!= k
    } // for k
} // for j
} // for i
} // for bords
hes << "Jacobian full" << endl;
for(int i = 0; i < pars.size() ; i++ ) {
    hes <<  setprecision(2)  << Jacobian[i] << '\t' ;
}
hes << endl;
hes << "Hessian full " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
    for(int k = 0; k < pars.size() ; k++ ) {
        hes <<  setprecision(2)  <<  Hessian[j][k] << '\t' ;
    }
    hes << endl;
}
occm.close();
gsl_matrix * IHessian = gsl_matrix_alloc(pars.size(),pars.size() );
const size_t N = pars.size();
gsl_permutation * p = gsl_permutation_alloc (N);
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
    Hvector = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
    gsl_matrix_set_row(IHessian , j , Hvector);
}
int oint =1;
int* spoint =&oint;
gsl_linalg_LU_decomp (IHessian,  p, spoint);
gsl_matrix * inverse=gsl_matrix_alloc(pars.size(),pars.size() );
gsl_linalg_LU_invert (IHessian, p, inverse);
hes << " inverse Hes " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
    gsl_matrix_get_row( Hvector ,inverse , j );
    for(int i = 0; i < pars.size() ; i++ ) {
        hes << gsl_vector_get(Hvector,i)  << '\t' ;
    }
    hes << endl;
}
double epsrel = 0 ; //.01;
gsl_matrix *covar;
covar = gsl_matrix_alloc( pars.size(), pars.size() );
gsl_multifit_covar(JM, epsrel, covar);
hes << "covar = JtJ inverse; error = sqrt(JtJ^-1) ,pg 446 gsl man" << endl;   /// why does this say inverse?
for( int i =0 ; i < covar->size1; i++) {
    for( int j =0 ; j < covar->size2; j++) {
        hes << sqrt(gsl_matrix_get(covar,i,j)) << '\t';
    }
    hes << endl;
}
hes << endl;
gsl_matrix * JMt = gsl_matrix_alloc( JM->size2, JM->size1);
gsl_matrix_transpose_memcpy(JMt,JM);   //(dest,source)
gsl_matrix * CV = gsl_matrix_alloc( JM->size2, JM->size2);
int oint2 =1;
int* spoint2 =&oint2;
gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
1.0, JMt, JM,
0.0, CV);    // dgemm simply multiples two matrices:  JMt*JM,.. namely J transpose time J, and stores in CV
gsl_linalg_LU_decomp (CV,  p, spoint2);
gsl_matrix * inverseCV=gsl_matrix_alloc( JM->size2, JM->size2);
gsl_linalg_LU_invert (CV, p, inverseCV);
hes << " inverse JtJ = Covar, manual " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(JM->size2 );
    gsl_matrix_get_row( Hvector ,inverseCV , j );
    for(int i = 0; i < JM->size2 ; i++ ) {
        hes << gsl_vector_get(Hvector,i)  << '\t' ;
    }
    hes << endl;
}
hes << endl;
hes << "  RMSE*(JtJ)^(-1) , parameter uncertainties , manual " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(JM->size2 );
    gsl_matrix_get_row( Hvector ,inverseCV , j );
    for(int i = 0; i < JM->size2 ; i++ ) {
        hes << jRMSE.getObj()* gsl_vector_get(Hvector,i) << '\t' ;
    }
    hes << endl;
}
hes  << " hes << jRMSE.getObj() " << jRMSE.getObj() << endl;
vector<int> indices(0);
for(int bords = 0; bords< AllData.size(); bords++) {  //bords represent vectors like siteVec.. each bord has all sequences sequence 		( for given bord, loop over all seqs)						// there are 12 bords, unless seqSitesmbot and 		the other 	snail vector are added.
for(int i=0; i < AllData[bords].size(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
indices.push_back(i);
}
}
hes.close();
cout << " Alldata.size " << AllData.size() << " AllBorders.size() "<<AllBorders.size()<< endl;
cout << " Alldata[0].size " << AllData[0].size() << " AllBorders[0].size() "<<AllBorders[0].size()<< endl;
cout << " Alldata[1].size " << AllData[1].size() << " AllBorders[1].size() "<<AllBorders[1].size()<< endl;
cout <<  " AllBorders[0] " << endl<<AllBorders[0]<< endl;
cout << "Alldata mat is mxn matrix while AllBorders mat is nxm mat." << endl;
for(int m = 0; m < seqSites.size(); m++ ) {
    tsites = seqSites[m];
    tsitesbot = seqSitesbot[m]  ;
    ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
    foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
    for( int i = 0; i < tsites.size() ; i++ ) {
        for (;;) {
            if ( j < tsites[i].start ) { foo << "b"; j++;  }
            else {//cout << " j = " << j << endl;
            if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
            if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
            if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
        }// else
    }
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {
        if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
        else {//cout << " j = " << j << endl;
        if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo.close();
j=0;
int jjj=10000;
for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
    if( i == 0 ) { fo << "&&" ; continue; }
    for( int ii = 0; ii < seqSitesm1[m].size(); ii++ ) {
        if ( i== seqSitesm1[m][ii].start ) {
            j=0;
            if ( abs(i-jjj) < 6 ) {
                while( j < 9 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                fo << "\\color{yellow}"<< ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++;}
                jjj=i; break;
            }//if abs
            else{
                if( seqSitesm1[m][ii].factorIdx == 0 ) {
                    while( j < 9 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                    fo << "\\color{blue}"<< ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++;}
                    jjj=i; break;
                }//if tsitesii
                if( seqSitesm1[m][ii].factorIdx == 1 ) {
                    while( j < 6 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                    fo << "\\color{green}"<<ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++; jjj=i;}
                    break;
                }  // although j > tsites[i++], hence they'll appear right after one another.
                if(seqSitesm1[m][ii].factorIdx == 2 ) {
                    while( j < 6 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                    fo << "\\color{red}"<<ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}";  i++; j++; jjj=i;}
                    break;
                }
            }// else abs()
        }// if
    }// for ii
    fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
    if(i% 90 == 0 ) {fo << "\\\\&&"; }
}//for i
j=0;
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
vector< double > concs = factorExprData.getCol( cell[2*m] );
double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< " " << setprecision(2) << predicted  << "&" << "cell "<< cell[2*m] << "&" ;
for( int i = 0; i < tsites.size() ; i++ ) {
    for (;;) {
        if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsites[i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {  // tsitesbot was changed from tsites below
    if ( j < tsitesbot[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
    else {//cout << " j = " << j << endl;
    if( j % 90 == 0 ) { fo << "\\\\&&"; }
    if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
    if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
    if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
}// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "d+" << "&" << "cell "<< cell[2*m]<< "&" ;
for( int i = 0; i < seqSitesf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesf2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "d+"  << "&" << "cell "<< cell[2*m+1]<< "&" ;
for( int i = 0; i < seqSitesbotf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesbotf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesbotf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesbotf2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesbotf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "m1"<< "d+" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1f2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] <<"m2"<< "d+" <<"&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2f2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\&&\\\\&&\\\\" << endl;
}// for m
fo.close();
}
void ExprPredictor::printFile3c( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) // came from printFile3b above 8/1/2025
{
    os.setf( ios::fixed );
    os.precision( 3 );
    os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    os << endl<<" txpEffects : "<< endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        os << par.txpEffects[i] << "\t";
    }
    os << endl<<  " coopertivity : " << endl;
    for(int i =0; i < nFactors(); i++){
        for(int j =0; j <= i; j++){
            if ( coopMat(i,j)){
                for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
                os <<  "\t" << par.theV[i][j][k] ;
            }
        }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
    }
}
os<< endl << " quenching : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <i; j++){
        if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
        for(int k =0; k < ExprPar::nbins ; k++){
            os <<  "\t" << par.theVr[i][j][k] ;
        }
    }
}
}
os << endl << "objective function : " << endl;
os << jRMSE.getObj() << endl;
ofstream fo( "oData/format.tex",ios::app );
int j=0;
vector< Site > tsites;
vector< Site > tsitesbot;
// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
double bottomofdborder = .5;
double topofdborder = .8;
int nrow = factorExprData.nRows();
int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
for(int m = 0; m < seqSites.size(); m++ ) {
    tsites = seqSites[m];
    tsitesbot = seqSitesbot[m]  ;
    ofstream foo( "format.txt",ios::app );
    foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
    for( int i = 0; i < tsites.size() ; i++ ) {
        for (;;) {
            if ( j < tsites[i].start ) { foo << "b"; j++;  }
            else {//cout << " j = " << j << endl;
            if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
            if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
            if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
        }// else
    }
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {
        if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
        else {//cout << " j = " << j << endl;
        if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo.close();
j=0;
int jjj=10000;
for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
    if( i == 0 ) { fo << "&&" ; continue; }
    for( int ii = 0; ii < seqSitesm1[m].size(); ii++ ) {
        if ( i== seqSitesm1[m][ii].start ) {
            j=0;
            if ( abs(i-jjj) < 6 ) {
                while( j < 9 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                fo << "\\color{yellow}"<< ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++;}
                jjj=i; break;
            }//if abs
            else{
                if( seqSitesm1[m][ii].factorIdx == 0 ) {
                    while( j < 9 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                    fo << "\\color{blue}"<< ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++;}
                    jjj=i; break;
                }//if tsitesii
                if( seqSitesm1[m][ii].factorIdx == 1 ) {
                    while( j < 6 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                    fo << "\\color{green}"<<ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}"; i++; j++; jjj=i;}
                    break;
                }  // although j > tsites[i++], hence they'll appear right after one another.
                if(seqSitesm1[m][ii].factorIdx == 2 ) {
                    while( j < 6 ) {if( i % 90 == 0 ) { fo << "\\\\&&"; }
                    fo << "\\color{red}"<<ALPHABET[ExprPredictor::seqsy[m][i]] <<"\\color{black}";  i++; j++; jjj=i;}
                    break;
                }
            }// else abs()
        }// if
    }// for ii
    fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
    if(i% 90 == 0 ) {fo << "\\\\&&"; }
}//for i
j=0;
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
vector< double > concs = factorExprData.getCol( cell[2*m] );
double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< " " << setprecision(2) << predicted  << "&" << "cell "<< cell[2*m] << "&" ;
for( int i = 0; i < tsites.size() ; i++ ) {
    for (;;) {
        if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsites[i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {  // tsitesbot was changed from tsites below
    if ( j < tsitesbot[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
    else {//cout << " j = " << j << endl;
    if( j % 90 == 0 ) { fo << "\\\\&&"; }
    if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
    if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
    if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
}// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "d+" << "&" << "cell "<< cell[2*m]<< "&" ;
for( int i = 0; i < seqSitesf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesf2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "d+"  << "&" << "cell "<< cell[2*m+1]<< "&" ;
for( int i = 0; i < seqSitesbotf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesbotf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesbotf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesbotf2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesbotf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "m1"<< "d+" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1f2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] <<"m2"<< "d+" <<"&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2f2[m][i].factorIdx == 1 ) { fo << "\\color{green}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\&&\\\\&&\\\\" << endl;
}// for m
fo.close();
}
void ExprPredictor::printFile3( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE )
{
    os.setf( ios::fixed );
    os.precision( 3 );
    os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    os << endl<<" txpEffects : "<< endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        os << par.txpEffects[i] << "\t";
    }
    os << endl<<  " coopertivity : " << endl;
    for(int i =0; i < nFactors(); i++){
        for(int j =0; j <= i; j++){
            if ( coopMat(i,j)){
                for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
                os <<  "\t" << par.theV[i][j][k] ;
            }
        }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
    }
}
os<< endl << " quenching : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <i; j++){
        if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
        for(int k =0; k < ExprPar::nbins ; k++){
            os <<  "\t" << par.theVr[i][j][k] ;
        }
    }
}
}
os << endl << "objective function : " << endl;
os << jRMSE.getObj() << endl;
ofstream fo( "oData/format.tex",ios::app );
int j=0;
vector< Site > tsites;
vector< Site > tsitesbot;
// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
double bottomofdborder = .5;
double topofdborder = .8;
int nrow = factorExprData.nRows();
int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
ofstream hes( "hes.txt" );
vector< double > pars;
par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators );
int pars_size = pars.size();
fix_pars.clear();
free_pars.clear();
for( int index = 0; index < pars_size; index++ ){
    if( indicator_bool[ index ] ){
        free_pars.push_back( pars[ index ]);
    }
    else{
        fix_pars.push_back( pars[ index ] );
    }
}
pars.clear();
pars = free_pars;
double step = .01;
bool a=1;
ofstream occm("occmat.txt");
cout << AllData.size() << endl;
cout << "AllData[0].size()"<< AllData[0].size() << endl;
cout << "AllData[1].size()"<< AllData[1].size() << endl;
cout << "AllBorders.size() " << AllBorders.size() << endl;
vector< double > initilizer( pars.size() );
Jacobian = initilizer;
vector< vector< double > > init2( pars.size(), initilizer );
Hessian = init2;
gsl_matrix *JM;  // number of [sequences * (number of bords)] by the number of parameters
int count = 0; // used to indicate data point within bords*seq
JM = gsl_matrix_alloc( AllData[0].size()* AllData.size(),pars.size() );  // number of sequences(*number of bords) by the number of parameters
for(int bords = 0; bords< AllData.size(); bords++) {  //bords represent vectors like siteVec.. each bord has all sequences sequence ( for given bord, loop over all seqs)
for(int i=0; i < AllData[bords].size(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
for (int jj = 0; jj < pars.size() ; jj++ ) {
    vector< double > parsjuh =pars;
    vector< double > parsjdh= pars;
    parsjuh[jj] = pars[jj] + step;
    parsjdh[jj] = pars[jj] - step;
    vector< double > concs = factorExprData.getCol( AllBorders[i][bords] );
    vector < double > all_pars;
    all_pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    for( int index = 0; index < indicator_bool.size(); index ++ ){
        if(  indicator_bool[ index ]  ){
            all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
        }
        else{
            all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
        }
    }
    ExprPar parjuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators , a);
    all_pars.clear();
    free_par_counter = 0;
    fix_par_counter = 0;
    for( int index = 0; index < indicator_bool.size(); index ++ ){
        if(  indicator_bool[ index ]  ){
            all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
        }
        else{
            all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
        }
    }
    ExprPar parjdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
    vector< double > fOcc(0);
    hes.setf( ios::fixed);
    hes.width(3);
    count = bords* nSeqs();
    ExprFunc* func = this->createExprFunc( parjuh );
    double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
    func = this->createExprFunc( parjdh );
    double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
    func = this->createExprFunc( par_model );
    if( jj==0){
        func->predictOcc( AllData[bords][i], 5, concs,  fOcc );
        ;
        occm << ExprPredictor::seqNmes[i]<< '\t'<<AllBorders[i][bords] <<'\t' << fOcc << endl;
    } //if i ==0
    double p = func->predictExpr(AllData[bords][i], 5, concs ); // AllData[bords][i] = seqSites[i]
    Jacobian[jj] += ( pjuh- pjdh )/(2*step);
    gsl_matrix_set( JM,count+ i, jj, ( pjuh- pjdh )/(2*step) );  // i is seqs, jj is parameters
    parsjuh.clear();
    parsjdh.clear();
}  // for jj
for (int j = 0; j < pars.size() ; j++ ) {
    for (int k = 0; k < pars.size() ; k++ ) {
        vector< double > parsjuh =pars;
        vector< double > parsjdh= pars;
        vector< double > parskuh =pars;  // this and the next aren't being used
        vector< double > parskdh= pars;
        vector< double > parsjkdh= pars;
        vector< double > parsjkuh =pars;
        vector< double > parsjukdh= pars;
        vector< double > parsjdkuh =pars;
        parsjuh[j] = pars[j] + step;    // f_j    first partial
        parsjdh[j] = pars[j] - step;    // f_j
        parskuh[k] = pars[k] + step;    // f_k
        parskdh[k] = pars[k] - step;     // f_k
        parsjkuh[j] = pars[j] + step;    //  f_kj   seconder order partial, first step in construction
        parsjkuh[k] = pars[k] + step;     // f_kj  second step in construction: f(x+h,y+h)
        parsjkdh[j] = pars[j] - step;    // f_kj  second order,..  f(k-h,j-h)
        parsjkdh[k] = pars[k] - step;
        parsjukdh[j] = pars[j] + step;     //  f_kj  : f(k-h, j+h)
        parsjukdh[k] = pars[k] - step;
        parsjdkuh[j] = pars[j] - step;     //  f_kj  : f(k+h, j-d )
        parsjdkuh[k] = pars[k] + step;
        vector< double > concs = factorExprData.getCol( AllBorders[i][bords] );
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parjuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parjdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parskuh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parskdh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parkdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
        if( j == k ) {
            ExprFunc* func = this->createExprFunc( parjuh );
            double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjdh );
            double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( par_model );
            double p = func->predictExpr(AllData[bords][i], 5, concs );
            Hessian[j][k] += ( pjuh -2*p + pjdh )/(step*step) ;
            parsjuh.clear();
            parsjdh.clear();
            parskuh.clear();
            parskdh.clear();
            parsjkuh.clear();
            parsjkdh.clear();
            parsjdkuh.clear();
            parsjukdh.clear();
        }
        if ( j !=k ) {
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjkuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjkdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjkdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjdkuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjdkuh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjukdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjukdh = ExprPar ( all_pars, coopMat, actIndicators, repIndicators, a );
            ExprFunc* func = this->createExprFunc( parjuh );
            double pjuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjdh );
            double pjdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( par_model );
            double p = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parkuh );
            double pkuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parkdh );
            double pkdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjkuh );
            double pjkuh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjkdh );
            double pjkdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjukdh );
            double pjukdh = func->predictExpr(AllData[bords][i], 5, concs );
            func = this->createExprFunc( parjdkuh );
            double pjdkuh = func->predictExpr(AllData[bords][i], 5, concs );
            Hessian[j][k] +=( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step);
            parsjuh.clear();
            parsjdh.clear();
            parskuh.clear();
            parskdh.clear();
            parsjkuh.clear();
            parsjkdh.clear();
        }  // if j!= k
    } // for k
} // for j
} // for i
} // for bords
hes << "Jacobian full" << endl;
for(int i = 0; i < pars.size() ; i++ ) {
    hes <<  setprecision(2)  << Jacobian[i] << '\t' ;
}
hes << endl;
hes << "Hessian full " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
    for(int k = 0; k < pars.size() ; k++ ) {
        hes <<  setprecision(2)  <<  Hessian[j][k] << '\t' ;
    }
    hes << endl;
}
occm.close();
gsl_matrix * IHessian = gsl_matrix_alloc(pars.size(),pars.size() );
const size_t N = pars.size();
gsl_permutation * p = gsl_permutation_alloc (N);
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
    Hvector = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
    gsl_matrix_set_row(IHessian , j , Hvector);
}
int oint =1;
int* spoint =&oint;
gsl_linalg_LU_decomp (IHessian,  p, spoint);
gsl_matrix * inverse=gsl_matrix_alloc(pars.size(),pars.size() );
gsl_linalg_LU_invert (IHessian, p, inverse);
hes << " inverse Hes " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
    gsl_matrix_get_row( Hvector ,inverse , j );
    for(int i = 0; i < pars.size() ; i++ ) {
        hes << gsl_vector_get(Hvector,i)  << '\t' ;
    }
    hes << endl;
}
double epsrel = 0 ; //.01;
gsl_matrix *covar;
covar = gsl_matrix_alloc( pars.size(), pars.size() );
gsl_multifit_covar(JM, epsrel, covar);
hes << "covar = JtJ inverse; error = sqrt(JtJ^-1) ,pg 446 gsl man" << endl;   /// why does this say inverse?
for( int i =0 ; i < covar->size1; i++) {
    for( int j =0 ; j < covar->size2; j++) {
        hes << sqrt(gsl_matrix_get(covar,i,j)) << '\t';
    }
    hes << endl;
}
hes << endl;
gsl_matrix * JMt = gsl_matrix_alloc( JM->size2, JM->size1);
gsl_matrix_transpose_memcpy(JMt,JM);   //(dest,source)
gsl_matrix * CV = gsl_matrix_alloc( JM->size2, JM->size2);
int oint2 =1;
int* spoint2 =&oint2;
gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
1.0, JMt, JM,
0.0, CV);    // dgemm simply multiples two matrices:  JMt*JM,.. namely J transpose time J, and stores in CV
gsl_linalg_LU_decomp (CV,  p, spoint2);
gsl_matrix * inverseCV=gsl_matrix_alloc( JM->size2, JM->size2);
gsl_linalg_LU_invert (CV, p, inverseCV);
hes << " inverse JtJ = Covar, manual " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(JM->size2 );
    gsl_matrix_get_row( Hvector ,inverseCV , j );
    for(int i = 0; i < JM->size2 ; i++ ) {
        hes << gsl_vector_get(Hvector,i)  << '\t' ;
    }
    hes << endl;
}
hes << endl;
hes << "  RMSE*(JtJ)^(-1) , parameter uncertainties , manual " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(JM->size2 );
    gsl_matrix_get_row( Hvector ,inverseCV , j );
    for(int i = 0; i < JM->size2 ; i++ ) {
        hes << jRMSE.getObj()* gsl_vector_get(Hvector,i) << '\t' ;
    }
    hes << endl;
}
hes  << " hes << jRMSE.getObj() " << jRMSE.getObj() << endl;
vector<int> indices(0);
for(int bords = 0; bords< AllData.size(); bords++) {  //bords represent vectors like siteVec.. each bord has all sequences sequence 		( for given bord, loop over all seqs)						// there are 12 bords, unless seqSitesmbot and 		the other 	snail vector are added.
for(int i=0; i < AllData[bords].size(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
indices.push_back(i);
}
}
hes.close();
cout << " Alldata.size " << AllData.size() << " AllBorders.size() "<<AllBorders.size()<< endl;
cout << " Alldata[0].size " << AllData[0].size() << " AllBorders[0].size() "<<AllBorders[0].size()<< endl;
cout << " Alldata[1].size " << AllData[1].size() << " AllBorders[1].size() "<<AllBorders[1].size()<< endl;
cout <<  " AllBorders[0] " << endl<<AllBorders[0]<< endl;
cout << "Alldata mat is mxn matrix while AllBorders mat is nxm mat." << endl;
for(int m = 0; m < seqSites.size(); m++ ) {
    tsites = seqSites[m];
    tsitesbot = seqSitesbot[m]  ;
    ofstream foo( "oData/format2.txt",ios::app );
    foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
    for( int i = 0; i < tsites.size() ; i++ ) {
        for (;;) {
            if ( j < tsites[i].start ) { foo << "b"; j++;  }
            else {//cout << " j = " << j << endl;
            if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
            if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
            if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
        }// else
    }
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {
        if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
        else {//cout << " j = " << j << endl;
        if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo.close();
j=0;
for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
    if( i == 0 ) { fo << "&&" ; continue; }
    fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
    if(i% 90 == 0 ) {fo << "\\\\&&"; }
}
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
vector< double > concs = factorExprData.getCol( cell[2*m] );
double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< " " << setprecision(2) << predicted  << "&" << "cell "<< cell[2*m] << "&" ;
for( int i = 0; i < tsites.size() ; i++ ) {
    for (;;) {
        if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {  // tsitesbot was changed from tsites below
    if ( j < tsitesbot[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
    else {//cout << " j = " << j << endl;
    if( j % 90 == 0 ) { fo << "\\\\&&"; }
    if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
    if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
    if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
}// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "d+" << "&" << "cell "<< cell[2*m]<< "&" ;
for( int i = 0; i < seqSitesf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "d+"  << "&" << "cell "<< cell[2*m+1]<< "&" ;
for( int i = 0; i < seqSitesbotf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesbotf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesbotf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesbotf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesbotf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "m1"<< "d+" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] <<"m2"<< "d+" <<"&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\&&\\\\&&\\\\" << endl;
}// for m
fo.close();
}
void ExprPredictor::printFile5( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE )
{
    ofstream hes( "hes5.txt" );
    vector< double > pars;
    par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators );
    int pars_size = pars.size();
    hes << " pars " << endl << pars << endl;
    fix_pars.clear();
    free_pars.clear();
    for( int index = 0; index < pars_size; index++ ){
        if( indicator_bool[ index ] ){
            free_pars.push_back( pars[ index ]);
        }
        else{
            fix_pars.push_back( pars[ index ] );
        }
    }
    //Hassan start:
    pars.clear();
    pars = free_pars;
    double step = .0000000001;
    bool a = 1;
    //spaceweights=0;
    //Hassan end
    //ofstream occm("occmat.txt");
    //cout << AllData.size() << endl;
    //cout << "AllData[1].size()"<< AllData[1].size() << endl;
    //cout << AllBorders.size() << endl;
    //cout << "Allb[1].size()"<< AllBorders[1].size() << endl;
    //cout << "Allb[1].size()"<< AllBorders<< endl;
    vector< double > initilizer( pars.size() );
    Jacobian = initilizer;
    vector< vector< double > > init2( pars.size(), initilizer );
    Hessian = init2;
    gsl_matrix *JM;  // number of [sequences * (number of bords)] by the number of parameters
    int count = 0; // used to indicate data point within bords*seq
    ///////////////////////////////////////////////////////
    JM = gsl_matrix_alloc( nSeqs(),pars.size() );  // number of sequences(*number of bords) by the number of
    ///////////////////
    for (int j = 0; j < pars.size() ; j++ ) {
        for (int k = 0; k < pars.size() ; k++ ) {
            vector< double > parsjuh =pars;
            vector< double > parsjdh= pars;
            vector< double > parskuh =pars;  // this and the next aren't being used
            vector< double > parskdh= pars;
            vector< double > parsjkdh= pars;
            vector< double > parsjkuh =pars;
            vector< double > parsjukdh= pars;
            vector< double > parsjdkuh =pars;
            parsjuh[j] = pars[j] + step;    // f_j    first partial
            parsjdh[j] = pars[j] - step;    // f_j
            parskuh[k] = pars[k] + step;    // f_k
            parskdh[k] = pars[k] - step;     // f_k
            parsjkuh[j] = pars[j] + step;    //  f_kj   seconder order partial, first step in construction
            parsjkuh[k] = pars[k] + step;     // f_kj  second step in construction: f(k+h,j+h)
            parsjkdh[j] = pars[j] - step;    // f_kj  second order,..  f(k-h,j-h)
            parsjkdh[k] = pars[k] - step;
            parsjukdh[j] = pars[j] + step;     //  f_kj  : f(k-h, j+h)
            parsjukdh[k] = pars[k] - step;
            parsjdkuh[j] = pars[j] - step;     //  f_kj  : f(k+h, j-h )
            parsjdkuh[k] = pars[k] + step;
            // k fixed j up
            vector < double > all_pars;
            all_pars.clear();
            int free_par_counter = 0;
            int fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
            //printPar( parjuh );
            ///// k fixed j down
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
            //  j fixed k up
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parskuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
            // j fixed k down
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parskdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parkdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
            ///// calculate fxx or fxy
            //cout.setf( ios::fixed );
            //   cout.precision( 3 );
            //     cout.width( 8 );
            //hes.setf( ios::fixed);
            //hes.width(3);
            if( j == k ) {
                //ExprFunc* func = this->createExprFunc( parjuh );
                double pjuh =  this->compRMSE2(parjuh); //func->predictExpr3(seqSitesm1[i], 5, concs );
                // func = this->createExprFunc( parjdh );
                double pjdh = this->compRMSE2(parjuh); // func->predictExpr3(seqSitesm1[i], 5, concs );
                // func = this->createExprFunc( par_model );
                double p =  this->compRMSE2(par_model);//func->predictExpr3(seqSitesm1[i], 5, concs );
                //hes << ( pjuh -2*p + pjdh )/(step*step) << '\t' ;
                Hessian[j][k] += ( pjuh -2*p + pjdh )/(step*step) ;
                parsjuh.clear();
                parsjdh.clear();
                parskuh.clear();
                parskdh.clear();
                parsjkuh.clear();
                parsjkdh.clear();
                parsjdkuh.clear();
                parsjukdh.clear();
            }
            if ( j !=k ) {
                //////////// both u or d
                all_pars.clear();
                free_par_counter = 0;
                fix_par_counter = 0;
                for( int index = 0; index < indicator_bool.size(); index ++ ){
                    if(  indicator_bool[ index ]  ){
                        all_pars.push_back( parsjkuh[ free_par_counter ++ ]  );
                    }
                    else{
                        all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                    }
                }
                ExprPar parjkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
                ///// step by 2
                all_pars.clear();
                free_par_counter = 0;
                fix_par_counter = 0;
                for( int index = 0; index < indicator_bool.size(); index ++ ){
                    if(  indicator_bool[ index ]  ){
                        all_pars.push_back( parsjkdh[ free_par_counter ++ ]  );
                    }
                    else{
                        all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                    }
                }
                ExprPar parjkdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
                ////////////  mixed
                all_pars.clear();
                free_par_counter = 0;
                fix_par_counter = 0;
                for( int index = 0; index < indicator_bool.size(); index ++ ){
                    if(  indicator_bool[ index ]  ){
                        all_pars.push_back( parsjdkuh[ free_par_counter ++ ]  );
                    }
                    else{
                        all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                    }
                }
                ExprPar parjdkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
                ///// step by 2
                all_pars.clear();
                free_par_counter = 0;
                fix_par_counter = 0;
                for( int index = 0; index < indicator_bool.size(); index ++ ){
                    if(  indicator_bool[ index ]  ){
                        all_pars.push_back( parsjukdh[ free_par_counter ++ ]  );
                    }
                    else{
                        all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                    }
                }
                ExprPar parjukdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators , a);
                //ExprFunc* func = this->createExprFunc( parjuh );
                double pjuh =this->compRMSE2(parjuh);// func->predictExpr3(seqSitesm1[i], 5, concs );
                // func = this->createExprFunc( parjdh );
                double pjdh = this->compRMSE2(parjdh);//func->predictExpr3(seqSitesm1[i], 5, concs );
                //func = this->createExprFunc( par_model );
                double p =this->compRMSE2(par_model);// func->predictExpr3(seqSitesm1[i], 5, concs );
                //  func = this->createExprFunc( parkuh );
                double pkuh = this->compRMSE2(parkuh);//func->predictExpr3(seqSitesm1[i], 5, concs );
                // func = this->createExprFunc( parkdh );
                double pkdh =this->compRMSE2(parkdh);// func->predictExpr3(seqSitesm1[i], 5, concs );
                //	func = this->createExprFunc( parjkuh );
                double pjkuh =this->compRMSE2(parjkuh);// func->predictExpr3(seqSitesm1[i], 5, concs );
                // func = this->createExprFunc( parjkdh );
                double pjkdh =this->compRMSE2(parjkdh);// func->predictExpr3(seqSitesm1[i], 5, concs );
                //func = this->createExprFunc( parjukdh );
                double pjukdh =this->compRMSE2(parjukdh);// func->predictExpr3(seqSitesm1[i], 5, concs );
                // func = this->createExprFunc( parjdkuh );
                double pjdkuh =this->compRMSE2(parjdkuh);// func->predictExpr3(seqSitesm1[i], 5, concs );
                printPar(parjkuh );
                printPar(parjkdh );
                printPar(parjdkuh );
                printPar(parjkdh );
                hes << pjkuh << " " << pjukdh << " " << pjdkuh << " " << pjkdh << endl;
                hes << ( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step) << '\t' << "delta " << 4*step*step << endl ;
                Hessian[j][k] +=( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step);
                parsjuh.clear();
                parsjdh.clear();
                parskuh.clear();
                parskdh.clear();
                parsjkuh.clear();
                parsjkdh.clear();
            }  // if j!= k
        } // for k
        //hes << endl;
    } // for j
    //hes<< endl  <<endl ;
    //} // for i
    hes << endl;
    hes << "Hessian full " << endl;
    for(int j = 0; j < pars.size() ; j++ ) {
        for(int k = 0; k < pars.size() ; k++ ) {
            hes <<  setprecision(6)  <<  Hessian[j][k] << '\t' ;
        }
        hes << endl;
    }
    gsl_matrix * IHessian = gsl_matrix_alloc(pars.size(),pars.size() );
    const size_t N = pars.size();
    gsl_permutation * p = gsl_permutation_alloc (N);
    for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
        Hvector = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
        gsl_vector_scale (Hvector, .5);
        gsl_matrix_set_row(IHessian , j , Hvector);
    }
    int oint =1;
    int* spoint =&oint;
    gsl_linalg_LU_decomp (IHessian,  p, spoint);
    gsl_matrix * inverse=gsl_matrix_alloc(pars.size(),pars.size() );
    gsl_linalg_LU_invert (IHessian, p, inverse);
    hes << " inverse .5*Hes, hence var of paramters " << endl;
    for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
        gsl_matrix_get_row( Hvector ,inverse , j );
        for(int i = 0; i < pars.size() ; i++ ) {
            hes << gsl_vector_get(Hvector,i)  << '\t' ;
        }
        hes << endl;
    }
    hes << " standard error of paramters (sqrt((.5*hess)^-1)) " << endl;
    for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
        gsl_matrix_get_row( Hvector ,inverse , j );
        for(int i = 0; i < pars.size() ; i++ ) {
            if(i==j ){
                hes <<" stderr par "<<i << " = " << sqrt(gsl_vector_get(Hvector,i)) << '\t';
            }
        }
        hes << endl;
    }
    hes  << " hes << jRMSE.getObj() " << jRMSE.getObj() << endl;
    gsl_matrix * Hessian2 = gsl_matrix_alloc(pars.size(),pars.size() );
    const size_t N2 = pars.size();
    gsl_permutation * p2 = gsl_permutation_alloc (N2);
    for(int j = 0; j < pars.size() ; j++ ) {
        gsl_vector * Hvector2 =gsl_vector_alloc(pars.size() );
        Hvector2 = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
        gsl_matrix_set_row(Hessian2 , j , Hvector2);
    }
    gsl_vector *eval = gsl_vector_alloc (pars.size());
    gsl_matrix *evec = gsl_matrix_alloc (pars.size(), pars.size());
    gsl_eigen_symmv_workspace * w =
    gsl_eigen_symmv_alloc (pars.size());
    gsl_eigen_symmv (Hessian2, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec,
    GSL_EIGEN_SORT_ABS_ASC);
    {
        int i;
        for (i = 0; i < pars.size(); i++)
        {
            double eval_i
            = gsl_vector_get (eval, i);
            gsl_vector_view evec_i
            = gsl_matrix_column (evec, i);
            printf ("eigenvalue = %g\n", eval_i);
            printf ("eigenvector = \n");
            gsl_vector_fprintf (stdout,
            &evec_i.vector, "%g");
        }
    }
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    hes.close();
}
void ExprPredictor::printFile4( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE )
{
    os.setf( ios::fixed );
    os.precision( 3 );
    os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    os << endl<<" txpEffects : "<< endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        os << par.txpEffects[i] << "\t";
    }
    os << endl<<  " coopertivity : " << endl;
    for(int i =0; i < nFactors(); i++){
        for(int j =0; j <= i; j++){
            if ( coopMat(i,j)){
                for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
                os <<  "\t" << par.theV[i][j][k] ;
            }
        }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
    }
}
os<< endl << " quenching : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <i; j++){
        if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
        for(int k =0; k < ExprPar::nbins ; k++){
            os <<  "\t" << par.theVr[i][j][k] ;
        }
    }
}
}
os << endl << "objective function : " << endl;
os << jRMSE.getObj() << endl;
os.close();
ofstream fo( "oData/format.tex",ios::app );
int j=0;
vector< Site > tsites;
vector< Site > tsitesbot;
// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
double bottomofdborder = .5;
double topofdborder = .8;
int nrow = factorExprData.nRows();
int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
ofstream hes( "hes4.txt" );
vector< double > pars;
par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators );
int pars_size = pars.size();
fix_pars.clear();
free_pars.clear();
for( int index = 0; index < pars_size; index++ ){
    if( indicator_bool[ index ] ){
        free_pars.push_back( pars[ index ]);
    }
    else{
        fix_pars.push_back( pars[ index ] );
    }
}
pars.clear();
pars = free_pars;
double step = .001;
bool a = 1;
vector< double > initilizer( pars.size() );
Jacobian = initilizer;
vector< vector< double > > init2( pars.size(), initilizer );
Hessian = init2;
gsl_matrix *JM;  // number of [sequences * (number of bords)] by the number of parameters
int count = 0; // used to indicate data point within bords*seq
JM = gsl_matrix_alloc( nSeqs(),pars.size() );  // number of sequences(*number of bords) by the number of
for(int i=0; i < nSeqs(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
for (int jj = 0; jj < pars.size() ; jj++ ) {
    vector< double > parsjuh =pars;
    vector< double > parsjdh= pars;
    parsjuh[jj] = pars[jj] + step;
    parsjdh[jj] = pars[jj] - step;
    vector < double > all_pars;
    all_pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    for( int index = 0; index < indicator_bool.size(); index ++ ){
        if(  indicator_bool[ index ]  ){
            all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
        }
        else{
            all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
        }
    }
    ExprPar parjuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators , a);
    all_pars.clear();
    free_par_counter = 0;
    fix_par_counter = 0;
    for( int index = 0; index < indicator_bool.size(); index ++ ){
        if(  indicator_bool[ index ]  ){
            all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
        }
        else{
            all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
        }
    }
    ExprPar parjdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
    vector< double > fOcc(0);
    hes.setf( ios::fixed);
    hes.width(3);
    count = nSeqs();
    double pjuh = this->compRMSE3(parjuh, i); // func->predictExpr3(seqSitesm1[i], 5, concs );
    double pjdh =  this->compRMSE3(parjdh, i); //func->predictExpr3(seqSitesm1[i], 5, concs );
    Jacobian[jj] += ( pjuh- pjdh )/(2.0*step);
    gsl_matrix_set( JM,i, jj, ( pjuh- pjdh )/(2.0*step) );  // i is seqs, jj is parameters
    parsjuh.clear();
    parsjdh.clear();
}  // for jj
for (int j = 0; j < pars.size() ; j++ ) {
    for (int k = 0; k < pars.size() ; k++ ) {
        vector< double > parsjuh =pars;
        vector< double > parsjdh= pars;
        vector< double > parskuh =pars;  // this and the next aren't being used
        vector< double > parskdh= pars;
        vector< double > parsjkdh= pars;
        vector< double > parsjkuh =pars;
        vector< double > parsjukdh= pars;
        vector< double > parsjdkuh =pars;
        parsjuh[j] = pars[j] + step;    // f_j    first partial
        parsjdh[j] = pars[j] - step;    // f_j
        parskuh[k] = pars[k] + step;    // f_k
        parskdh[k] = pars[k] - step;     // f_k
        parsjkuh[j] = pars[j] + step;    //  f_kj   seconder order partial, first step in construction
        parsjkuh[k] = pars[k] + step;     // f_kj  second step in construction: f(k+h,j+h)
        parsjkdh[j] = pars[j] - step;    // f_kj  second order,..  f(k-h,j-h)
        parsjkdh[k] = pars[k] - step;
        parsjukdh[j] = pars[j] + step;     //  f_kj  : f(k-h, j+h)
        parsjukdh[k] = pars[k] - step;
        parsjdkuh[j] = pars[j] - step;     //  f_kj  : f(k+h, j-h )
        parsjdkuh[k] = pars[k] + step;
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parjuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parjdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parskuh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
        all_pars.clear();
        free_par_counter = 0;
        fix_par_counter = 0;
        for( int index = 0; index < indicator_bool.size(); index ++ ){
            if(  indicator_bool[ index ]  ){
                all_pars.push_back( parskdh[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar parkdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
        if( j == k ) {
            double pjuh =  this->compRMSE3(parjuh, i); //func->predictExpr3(seqSitesm1[i], 5, concs );
            double pjdh = this->compRMSE3(parjdh, i); // func->predictExpr3(seqSitesm1[i], 5, concs );
            double p =  this->compRMSE3(par_model, i);//func->predictExpr3(seqSitesm1[i], 5, concs );
            Hessian[j][k] += ( pjuh -2*p + pjdh )/(4*step*step) ;
            parsjuh.clear();
            parsjdh.clear();
            parskuh.clear();
            parskdh.clear();
            parsjkuh.clear();
            parsjkdh.clear();
            parsjdkuh.clear();
            parsjukdh.clear();
        }
        if ( j !=k ) {
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjkuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjkdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjkdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjdkuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjdkuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjukdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjukdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators , a);
            double pjuh =this->compRMSE3(parjuh, i);// func->predictExpr3(seqSitesm1[i], 5, concs );
            double pjdh = this->compRMSE3(parjdh, i);//func->predictExpr3(seqSitesm1[i], 5, concs );
            double p =this->compRMSE3(par_model, i);// func->predictExpr3(seqSitesm1[i], 5, concs );
            double pkuh = this->compRMSE3(parkuh,i);//func->predictExpr3(seqSitesm1[i], 5, concs );
            double pkdh =this->compRMSE3(parkdh, i);// func->predictExpr3(seqSitesm1[i], 5, concs );
            double pjkuh =this->compRMSE3(parjkuh, i);// func->predictExpr3(seqSitesm1[i], 5, concs );
            double pjkdh =this->compRMSE3(parjkdh, i);// func->predictExpr3(seqSitesm1[i], 5, concs );
            double pjukdh =this->compRMSE3(parjukdh, i);// func->predictExpr3(seqSitesm1[i], 5, concs );
            double pjdkuh =this->compRMSE3(parjdkuh, i);// func->predictExpr3(seqSitesm1[i], 5, concs );
            hes << pjkuh << " " << pjukdh << " " << pjdkuh << " " << pjkdh << endl;
            hes << ( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step) << '\t' << "delta " << 4*step*step << endl ;
            Hessian[j][k] +=( pjkuh -pjukdh - pjdkuh + pjkdh )/(4*step*step);
            parsjuh.clear();
            parsjdh.clear();
            parskuh.clear();
            parskdh.clear();
            parsjkuh.clear();
            parsjkdh.clear();
        }  // if j!= k
    } // for k
} // for j
} // for i
hes << "Jacobian full" << endl;
for(int i = 0; i < pars.size() ; i++ ) {
    hes <<  setprecision(2)  << Jacobian[i] << '\t' ;
}
hes << endl;
hes << "Hessian full " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
    for(int k = 0; k < pars.size() ; k++ ) {
        hes <<  setprecision(6)  <<  Hessian[j][k] << '\t' ;
    }
    hes << endl;
}
gsl_matrix * IHessian = gsl_matrix_alloc(pars.size(),pars.size() );
const size_t N = pars.size();
gsl_permutation * p = gsl_permutation_alloc (N);
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
    Hvector = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
    cout <<" gslvector before scale " << gsl_vector_get(Hvector,j) << endl ;
    gsl_vector_scale (Hvector, .5);
    cout <<" gslvector after scale " << gsl_vector_get(Hvector,j) << endl ;
    gsl_matrix_set_row(IHessian , j , Hvector);
}
int oint =1;
int* spoint =&oint;
gsl_linalg_LU_decomp (IHessian,  p, spoint);
gsl_matrix * inverse=gsl_matrix_alloc(pars.size(),pars.size() );
gsl_linalg_LU_invert (IHessian, p, inverse);
hes << " inverse .5*Hes, hence var of paramters " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
    gsl_matrix_get_row( Hvector ,inverse , j );
    for(int i = 0; i < pars.size() ; i++ ) {
        hes << gsl_vector_get(Hvector,i)  << '\t' ;
    }
    hes << endl;
}
hes << " standard error of paramters (sqrt((.5*hess)^-1)) " << endl;
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(pars.size() );
    gsl_matrix_get_row( Hvector ,inverse , j );
    for(int i = 0; i < pars.size() ; i++ ) {
        if(i==j ){
            hes <<" stderr par "<<i << " = " << sqrt(gsl_vector_get(Hvector,i)) << '\t';
        }
    }
    hes << endl;
}
double epsrel =.001;
gsl_matrix *covar;
covar = gsl_matrix_alloc( pars.size(), pars.size() );
gsl_multifit_covar(JM, epsrel, covar);
hes << "covar = JtJ inverse, sqrt((JtJ)^-1) " << endl;   /// why does this say inverse?
for( int i =0 ; i < covar->size1; i++) {
    for( int j =0 ; j < covar->size2; j++) {
        if(i==j ){
            hes <<" stderr par "<<i << " = " << sqrt(gsl_matrix_get(covar,i,j)) << '\t';
        }
    }
    hes << endl;
}
hes << endl;
gsl_matrix * JMt = gsl_matrix_alloc( JM->size2, JM->size1);
gsl_matrix_transpose_memcpy(JMt,JM);   //(dest,source)
gsl_matrix * CV = gsl_matrix_alloc( JM->size2, JM->size2);
int oint2 =1;
int* spoint2 =&oint2;
gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
1.0, JMt, JM,
0.0, CV);    // dgemm simply multiples two matrices:  JMt*JM,.. namely J transpose time J, and stores in CV
gsl_linalg_LU_decomp (CV,  p, spoint2);
gsl_matrix * inverseCV=gsl_matrix_alloc( JM->size2, JM->size2);
gsl_linalg_LU_invert (CV, p, inverseCV);
hes << " Jacobian " << endl;
for(int j = 0; j < JM->size1 ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(JM->size2 );
    gsl_matrix_get_row( Hvector ,JM , j );
    for(int i = 0; i < JM->size2 ; i++ ) {
        hes << gsl_vector_get(Hvector,i)  << '\t' ;
    }
    hes << endl;
}
hes << endl;
hes << " inverse JtJ = Covar, manual " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(JM->size2 );
    gsl_matrix_get_row( Hvector ,inverseCV , j );
    for(int i = 0; i < JM->size2 ; i++ ) {
        hes << gsl_vector_get(Hvector,i)  << '\t' ;
    }
    hes << endl;
}
hes << endl;
hes << "  MSE*(JtJ)^(-1) , parameter uncertainties , manual " << endl;
hes << "  MSE=sum_i ( y_i - ym(x_i) )^2 / (n-p) " << endl;
hes << "  where n= nSeqs()*nConds(), and p is the number of fit parameters " << endl;
hes << "  it's possible n should not include all nConds(), since they are/seem dependent " << endl;
for(int j = 0; j < JM->size2 ; j++ ) {
    gsl_vector * Hvector =gsl_vector_alloc(JM->size2 );
    gsl_matrix_get_row( Hvector ,inverseCV , j );
    for(int i = 0; i < JM->size2 ; i++ ) {
        hes << jRMSE.getObj() /( nSeqs()*nConds() )* gsl_vector_get(Hvector,i) << '\t' ;
    }
    hes << endl;
}
hes  << " chi^2jRMSE.getObj() " << jRMSE.getObj() << endl;
hes  << " hes << nSeqs()*nConds()" << nSeqs()*nConds() << endl;
hes  << " hes << chi^2/nSeqs()*nConds()" << jRMSE.getObj()/(nSeqs()*nConds()) << endl;
hes  << " hes << jRMSE.getObj() " << jRMSE.getObj() << endl;
gsl_matrix * Hessian2 = gsl_matrix_alloc(pars.size(),pars.size() );
const size_t N2 = pars.size();
gsl_permutation * p2 = gsl_permutation_alloc (N2);
for(int j = 0; j < pars.size() ; j++ ) {
    gsl_vector * Hvector2 =gsl_vector_alloc(pars.size() );
    Hvector2 = vector2gsl(Hessian[j]);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
    gsl_matrix_set_row(Hessian2 , j , Hvector2);
}
gsl_vector *eval = gsl_vector_alloc (pars.size());
gsl_matrix *evec = gsl_matrix_alloc (pars.size(), pars.size());
gsl_eigen_symmv_workspace * w =
gsl_eigen_symmv_alloc (pars.size());
gsl_eigen_symmv (Hessian2, eval, evec, w);
gsl_eigen_symmv_free (w);
gsl_eigen_symmv_sort (eval, evec,
GSL_EIGEN_SORT_ABS_ASC);
{
    int i;
    for (i = 0; i < pars.size(); i++)
    {
        double eval_i
        = gsl_vector_get (eval, i);
        gsl_vector_view evec_i
        = gsl_matrix_column (evec, i);
        printf ("eigenvalue = %g\n", eval_i);
        printf ("eigenvector = \n");
        gsl_vector_fprintf (stdout,
        &evec_i.vector, "%g");
    }
}
gsl_vector_free (eval);
gsl_matrix_free (evec);
hes.close();
}
void ExprPredictor::printFile2( ofstream& os, const ExprPar& par , ExprPredictor& jRMSE ) const
{
    os.setf( ios::fixed );
    os.precision( 3 );
    os << " maxBindingWts : " << endl;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            os << par.maxBindingWts[i] << "\t";
        }
    }
    os << endl<<" txpEffects : "<< endl;
    for ( int i = 0; i < nFactors(); i++ ) {
        os << par.txpEffects[i] << "\t";
    }
    os << endl<<  " coopertivity : " << endl;
    for(int i =0; i < nFactors(); i++){
        for(int j =0; j <= i; j++){
            if ( coopMat(i,j)){
                for(int k =0; k < par.theV[i][j].size(); k++) { //  i don't think this should be k less then vij.size...
                os <<  "\t" << par.theV[i][j][k] ;
            }
        }	//<< "\t" << "coopBin" <<"\t" << k << "\t"<<"factor" << i << j <<endl;
    }
}
os<< endl << " quenching : " << endl;
for(int i =0; i < nFactors(); i++){
    for(int j =0; j <i; j++){
        if(repIndicators[i])  { cout << endl;          // why is this cout here 9 13 11
        for(int k =0; k < ExprPar::nbins ; k++){
            os <<  "\t" << par.theVr[i][j][k] ;
        }
    }
}
}
os << endl << "objective function : " << endl;
os << jRMSE.getObj() << endl;
ofstream fo( "oData/format.tex",ios::app );
int j=0;
vector< Site > tsites;
vector< Site > tsitesbot;
// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
double bottomofdborder = .5;
double topofdborder = .8;
int nrow = factorExprData.nRows();
int ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
for(int m = 0; m < seqSites.size(); m++ ) {
    tsites = seqSites[m];
    tsitesbot = seqSitesbot[m]  ;
    ofstream foo( "/home/jacobc/Desktop/Desktop1026/clustalw/format2.txt",ios::app );
    foo<< ">" << ExprPredictor::seqNmes[m]  <<endl ;
    for( int i = 0; i < tsites.size() ; i++ ) {
        for (;;) {
            if ( j < tsites[i].start ) { foo << "b"; j++;  }
            else {//cout << " j = " << j << endl;
            if( tsites[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
            if( tsites[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
            if( tsites[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
        }// else
    }
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo<< ">" << ExprPredictor::seqNmes[m] <<"bot" <<endl ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {
        if ( j < tsitesbot[i].start ) { foo << "b"; j++;  }
        else {//cout << " j = " << j << endl;
        if( tsitesbot[i].factorIdx == 0 ) { foo << "d"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsitesbot[i].factorIdx == 1 ) { foo << "t"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsitesbot[i].factorIdx == 2 ) { foo << "s"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {
    foo << "b"; j++;
}
foo << endl;
foo.close();
j=0;
for(int i = 0; i < ExprPredictor::seqsy[m].size(); i++ ) {
    if( i == 0 ) { fo << "&&" ; continue; }
    fo << ALPHABET[ExprPredictor::seqsy[m][i]] ;
    if(i% 90 == 0 ) {fo << "\\\\&&"; }
}
fo << "\\\\&&\\\\";
ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
vector< double > concs = factorExprData.getCol( cell[2*m] );
double predicted = func->predictExpr(tsites, 5, concs );
fo << ExprPredictor::seqNmes[m]<< " " << predicted  << "&" << "pcell "<< jRMSE.cell[2*m] << "&" ;
for( int i = 0; i < tsites.size() ; i++ ) {
    for (;;) {
        if ( j < tsites[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( tsites[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( tsites[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( tsites[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "&" << "cell "<< cell[2*m+1] << "&" ;
for( int i = 0; i < tsitesbot.size() ; i++ ) {
    for (;;) {  // tsitesbot was changed from tsites below
    if ( j < tsitesbot[i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
    else {//cout << " j = " << j << endl;
    if( j % 90 == 0 ) { fo << "\\\\&&"; }
    if( tsitesbot[i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
    if( tsitesbot[i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
    if( tsitesbot[i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
}// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<<"m1" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "m2" << "&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "d+" << "&" << "cell "<< cell[2*m]<< "&" ;
for( int i = 0; i < seqSitesf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m]<< "d+"  << "&" << "cell "<< cell[2*m+1]<< "&" ;
for( int i = 0; i < seqSitesbotf2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesbotf2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesbotf2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesbotf2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesbotf2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] << "m1"<< "d+" << "&" << "cell "<< gsl_vector_get(transition_Indicests,2) << "&" ;
for( int i = 0; i < seqSitesm1f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm1f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm1f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm1f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm1f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\";
j=0;
fo << ExprPredictor::seqNmes[m] <<"m2"<< "d+" <<"&" << "cell "<< gsl_vector_get(transition_Indicesbs,2)<< "&" ;
for( int i = 0; i < seqSitesm2f2[m].size() ; i++ ) {
    for (;;) {
        if ( j < seqSitesm2f2[m][i].start ) { fo << "e"; j++; if( j % 90 == 0 ) { fo << "\\\\&&"; }}
        else {//cout << " j = " << j << endl;
        if( j % 90 == 0 ) { fo << "\\\\&&"; }
        if( seqSitesm2f2[m][i].factorIdx == 0 ) { fo << "\\color{blue}{d}\\color{black}"; j++;  break; }  // if two sites start at same position we shouldn't increment j
        if( seqSitesm2f2[m][i].factorIdx == 1 ) { fo << "\\color{blue}{t}\\color{black}"; j++ ; break; }  // although j > tsites[i++], hence they'll appear right after one another.
        if( seqSitesm2f2[m][i].factorIdx == 2 ) { fo << "\\color{red}{s}\\color{black}"; j++ ; break; }
    }// else
}
}
while( j < ExprPredictor::seqsy[m].size() ) {if( j % 90 == 0 ) { fo << "\\\\&&"; }
fo << "e"; j++;
}
fo << "\\\\&&\\\\&&\\\\" << endl;
}// for m
fo.close();
}
int ExprPredictor::train()
{
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;	// create rng type
    rng = gsl_rng_alloc( T );
    gsl_rng_set( rng, 12345 );		// set the seed equal to simulTime(0)
    ExprPar par_default( nFactors() );
    train( par_default, rng );
    return 0;
}
int ExprPredictor::predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs )
{
    targetExprs.clear();
    ExprFunc* func = this->createExprFunc( par_model );	// changed from createExprFunc2 on 92511
    for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr3( targetSites, targetSeqLength, concs );
        targetExprs.push_back( predicted );
    }
    return 0;
}
ModelType ExprPredictor::modelOption = LOGISTIC;
int ExprPredictor::estBindingOption = 1;    // 1. estimate binding parameters; 0. not estimate binding parameters
ObjType ExprPredictor::objOption = SSE;
double ExprPredictor::exprSimCrossCorr( const vector< double >& x, const vector< double >& y )
{
    vector< int > shifts;
    for ( int s = -maxShift; s <= maxShift; s++ ) {
        shifts.push_back( s );
    }
    vector< double > cov;
    vector< double > corr;
    cross_corr( x, y, shifts, cov, corr );
    double result = 0, weightSum = 0;
    result = *max_element( corr.begin(), corr.end() );
    return result;
}
int ExprPredictor::maxShift = 5;
double ExprPredictor::shiftPenalty = 0.8;
int ExprPredictor::nAlternations = 4;
int ExprPredictor::nRandStarts = 0;
double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
int ExprPredictor::nSimplexIters = 3;
int ExprPredictor::nGradientIters = 6;
int ExprPredictor::randSamplePar( const gsl_rng* rng, ExprPar& par ) const
{
    int counter= 0 ;
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if( indicator_bool[counter] ) {
                double rand_weight = exp( gsl_ran_flat( rng, log( ExprPar::min_weight ), log( ExprPar::max_weight ) ) );
                par.maxBindingWts[i] = rand_weight;
                counter++;
            }
            else{ counter++; }
        }
    }
    if(modelOption == BINS){
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                if ( coopMat( i, j ) ) {
                    for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
                    if( indicator_bool[counter] ) {
                        double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interaction ), log( ExprPar::max_interaction ) ) );
                        par.theV[ i][ j][k] = rand_interaction;
                        counter++;
                    }       //for
                    else{counter++;}
                } //if indicator bool
            }  //if construct
        }
    }
    for ( int i = 0; i < nFactors(); i++ ) {
        for ( int j = 0; j < i; j++ ) {
            if (repIndicators[i] ) {
                for(int k=0; k< ExprPar::nbins; k++){  // i.numbins, j.numbins..
                if( indicator_bool[counter] ) {
                    double rand_interaction = exp( gsl_ran_flat( rng, log( ExprPar::min_interactionr ), log( ExprPar::max_interactionr ) ) );
                    par.theVr[ i][ j][k] = rand_interaction;
                    counter++;
                }// for k
                else{counter++;}
            }// if indicator
        }  //if construct
    }
}
}
for ( int i = 0; i < nFactors(); i++ ) {
    if( indicator_bool[counter] ) {
        double rand_effect =  gsl_ran_flat( rng,  ExprPar::min_effect_Thermo , ExprPar::max_effect_Thermo  );
        par.txpEffects[i] = rand_effect;
        counter++;
    }
    else{counter++;}
}
return 0;
}
bool ExprPredictor::testPar( const ExprPar& par ) const
{
    if ( estBindingOption ) {
        for ( int i = 0; i < nFactors(); i++ ) {
            if ( par.maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) || par.maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) )
            return false;
        }
    }
    if( modelOption == BINS) {
        for ( int i = 0; i < nFactors(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                if(coopMat(i,j)){
                    for(int k =0; k < ExprPar::nbins; k++){
                        if ( par.theV[i][j][k] < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) || par.theV[i][j][k] > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) ) {
                            return false; }
                        }
                    }
                }
            }
            for ( int i = 0; i < nFactors(); i++ ) {
                for ( int j = 0; j <= i; j++ ) {
                    if(repIndicators[i]){
                        for(int k =0; k < ExprPar::nbins; k++){
                            if ( par.theVr[i][j][k] < ExprPar::min_interactionr * ( 1.0 + ExprPar::delta ) || par.theVr[i][j][k] > ExprPar::max_interactionr * ( 1.0 - ExprPar::delta ) ) {
                                return false; }
                            }
                        }
                    }
                }
            }
            for ( int i = 0; i < nFactors(); i++ ) {
                if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
            }
            return true;
        }
        void ExprPredictor::printPar2(  )
        {
            cout <<"dltwrep" << endl;
        }
        void ExprPredictor::printPar( const ExprPar& par ) const
        {
            cout.setf( ios::fixed );
            cout.precision( 3 );
            //     cout.width( 8 );
            // print binding weights
            if ( estBindingOption ) {
                for ( int i = 0; i < nFactors(); i++ ) {
                    if(i == 0){
                        cout << par.maxBindingWts[i] << "\t";
                    }
                    if(i==1){
                        cout << par.maxBindingWts[i] << "\t";
                    }
                    if(i==2){
                        cout << par.maxBindingWts[i] << "\t";
                    }
                }
            }
            //if(modelOption== BINS){
                for(int i =0; i < nFactors(); i++){
                    for(int j =0; j <= i; j++){
                        if(coopMat(i,j))  { cout << endl;
                        for(int k =0; k < ExprPar::nbins ; k++){   // was an erro 525 had k < theV(ij)
                        if( i == 1 && j ==0 ) {
                            cout << par.theV[i][j][k]  <<"\t";}
                            //count = count + 1 ;//ouut << endl;
                            if( i == 2 && j ==0 ) {
                                cout  << par.theV[i][j][k]  <<"\t";}
                                //count = count + 1 ;//ouut << endl;
                                if( i== 2 && j == 1) {
                                    cout  << par.theV[i][j][k]  <<"\t";}
                                    //count = count + 1 ;//ouut << endl;
                                }
                            }
                        }
                    }
                    for(int i =0; i < nFactors(); i++){
                        for(int j =0; j <=i; j++){
                            if(repIndicators[i])  { cout << endl;
                            for(int k =0; k < ExprPar::nbins ; k++){   // was an erro 525 had k < theV(ij)
                            if( i == 1 && j ==0 ) {
                                cout << par.theVr[i][j][k]  <<"\t";}
                                //count = count + 1 ;//ouut << endl;
                                if( i == 2 && j == 0) {
                                    cout << par.theVr[i][j][k]  <<"\t";}
                                    //count = count + 1 ;//ouut << endl;
                                    if( i== 2 && j == 1) {
                                        cout  << par.theVr[i][j][k]  <<"\t";}
                                        //count = count + 1 ;//ouut << endl;
                                    }
                                }
                            }
                        }
                        //}
                        // print the transcriptional effects
                        for ( int i = 0; i < nFactors(); i++ ) {
                            cout << par.txpEffects[i] << "\t";
                        }
                        cout << endl;
                    }
                    ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
                    {
                        return new ExprFunc( motifs, intFunc, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, par,  binBord, binBordr );
                    }
                    ExprFunc* ExprPredictor::createExprFunc2(  ExprPar& par )  // the heading of an accessor method is followed by const, to tell the compiler that such a method may not modify any of the instance variables.,  a declaration ,starting with const, (e.g. const ExprPar& Par) of a named variable or parameter (e.g. par) indicates that this variable can not be modified..  by declaring a parametere whose type is type& we have the type object refer directly to the object instead of making a copy of the object being passed,  this is called a reference parameter..
                    {
                        return new ExprFunc( motifs, intFunc, actIndicators, maxContact, repIndicators, repressionMat, repressionDistThr, par, binBord , binBordr);
                    }// motifs, intFunc.. are all private members of ExprPredictor... whereas par is passed as a reference parameter..
                    double ExprPredictor::compRMSE( const ExprPar& par ) const
                    {
                        ExprFunc* func = createExprFunc( par );
                        double rss = 0;
                        double squaredErr = 0;
                        for ( int i = 0; i < nSeqs(); i++ ) {
                            vector< double > predictedExprs;
                            vector< double > observedExprs;
                            for ( int j = 0; j < nConds(); j++ ) {
                                vector< double > concs = factorExprData.getCol( j );
                                double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
                                predictedExprs.push_back( predicted );
                                double observed = exprData( i, j );
                                observedExprs.push_back( observed );
                            }
                            double beta;
                            squaredErr += least_square( predictedExprs, observedExprs, beta );
                            for ( int i = 0; i < nConds(); i++ ) {  // 10/04
                            rss += ( predictedExprs[i] -  observedExprs[i] ) * ( predictedExprs[i] -  observedExprs[i] ); //beta*x is p, where x is now acting as a above, and beta is x above.
                        }
                    }
                    double rmse =  rss/ ( nSeqs() * nConds() ) ;
                    return rmse;
                }
                double ExprPredictor::compf( const ExprPar& par ,  vector<double>&ff)  //gsl_vector * f )
                {
                    spaceweights=0;
                    ExprFunc* func = createExprFunc( par );
                    double rss = 0;
                    double genespaceweight=0;
                    for ( int i = 0; i < nSeqs(); i++ ) {
                        vector< double > predictedExprs;
                        vector< double > observedExprs;
                        for ( int j = 0; j < nConds(); j++ ) {
                            vector< double > concs = factorExprData.getCol( j );
                            double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
                            predictedExprs.push_back( predicted );
                            double observed = exprData( i, j );
                            observedExprs.push_back( observed );
                        }
                        vector< double > concsweight(0);
                        for ( int j = 1; j < nConds()-1; j++ ) {
                            concsweight.push_back(  (exprData(i,j )-exprData(i,j+1) )/ 2  + (exprData(i,j) - exprData(i,j-1) )/2 );
                            genespaceweight=genespaceweight+concsweight[j];
                        }
                        for ( int j = 0; j < nConds(); j++ ) {  // 10/04
                        double obser = observedExprs[j];
                        double predic =predictedExprs[j];
                        if(observedExprs[j]==0){ obser =.01;}
                        if(observedExprs[j]==1){ obser =.99;}
                        rss +=( predic -  obser) * ( predic -  obser ); // / (predic*(1-predic ));
                        ff.push_back( (predic-obser)/ sqrt(obser));
                    }
                }
                for ( int i = 0; i < nSeqs(); i++ ) {
                    for ( int j = 0; j < nConds(); j++ ) {
                    }
                }
                double rmse = rss;
                return rmse; // /(nSeqs()*nConds() );
            }
            double ExprPredictor::compRMSE2( const ExprPar& par )
            {
                spaceweights=0;
                ExprFunc* func = createExprFunc( par );
                double rss = 0;
                double genespaceweight=0;
                for ( int i = 0; i < nSeqs(); i++ ) {
                    vector< double > predictedExprs;
                    vector< double > observedExprs;
                    for ( int j = 0; j < nConds(); j++ ) {
                        vector< double > concs = factorExprData.getCol( j );
                        double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
                        predictedExprs.push_back( predicted );
                        double observed = exprData( i, j );
                        observedExprs.push_back( observed );
                    }
                    vector< double > concsweight(0);
                    for ( int j = 1; j < nConds()-1; j++ ) {
                        concsweight.push_back(  (exprData(i,j )-exprData(i,j+1) )/ 2  + (exprData(i,j) - exprData(i,j-1) )/2 );
                        genespaceweight=genespaceweight+concsweight[j];
                    }
                    for ( int j = 0; j < nConds(); j++ ) {  // 10/04
                    double obser = observedExprs[j];
                    double predic =predictedExprs[j];
                    rss +=( predic -  obser) * ( predic -  obser ); // / (predic*(1-predic ));
                }
            }
            spaceweights=genespaceweight;
            double rmse =  rss;   //  / ( nSeqs() * nConds() )  ;
            return rmse; // /(nSeqs()*nConds() );
        }
        double ExprPredictor::compRMSE3( const ExprPar& par , int i)
        {
            ExprFunc* func = createExprFunc( par );
            double rss = 0;
            double squaredErr = 0;
            vector< double > predictedExprs;
            vector< double > observedExprs;
            for ( int j = 0; j < nConds(); j++ ) {
                vector< double > concs = factorExprData.getCol( j );
                double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
                predictedExprs.push_back( predicted );
                double observed = exprData( i, j );
                observedExprs.push_back( observed );
            }
            vector< double > concsweight(0);
            double genespaceweight=1; //1.0/nSeqs();  //setting to 1.0/nSeqs() is a pseudocount.
            for ( int j = 1; j < nConds()-1; j++ ) {
                concsweight.push_back(  (exprData(i,j )-exprData(i,j+1) )/ 2  + (exprData(i,j) - exprData(i,j-1) )/2 );
            }
            for ( int j = 0; j < nConds(); j++ ) {  // 10/04
            double obser = observedExprs[j];
            double predic =predictedExprs[j];
            rss +=( predic -  obser) * ( predic -  obser ); // / (predic*(1-predic ));
        }
        double rmse = rss; ///rss/ (  nConds() )  ;
        return   rmse;//  /( nSeqs() * nConds() ) ;
    }
    double ExprPredictor::compRMSE4( const ExprPar& par , int i, int j)
    {
        ExprFunc* func = createExprFunc( par );
        double rss = 0;
        double squaredErr = 0;
        vector< double > predictedExprs;
        vector< double > observedExprs;
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr3( seqSitesm1[ i ], seqLengths[i], concs );
        predictedExprs.push_back( predicted );
        double observed = exprData( i, j );
        observedExprs.push_back( observed );
        double obser=observed;
        if(observed==0){ obser =.01;}
        if(observed==1){ obser =.99;}
        return    predicted/sqrt(obser) ; /// (obser*(1-obser ));//rmse;
    }
    double ExprPredictor::compAvgCorr( const ExprPar& par ) const
    {
        ExprFunc* func = createExprFunc( par );
        vector< double > fOcc;
        vector< double > corrs;
        for ( int i = 0; i < seqSitesb.size(); i++ ) {   //seqs is a vector<vector<Sequence>> a Sequence has a DNA sequence and DNA datalike length
        vector< double > predicted;
        for ( int j = 0; j < seqSitesb[ i ].size(); j++ ) {
            vector< double > concs = factorExprData.getCol( 0 );  //103
            func->predictOcc( seqSitesb[ i ][ j ], i, concs,fOcc );
            predicted.push_back( fOcc[i] );// fOcc[i] represents the Experiment number, so if data for two experiments i will be 0 and 1.
        }
        corrs.push_back( correlation( predicted, bindingData[ i ] ) );
    }
    return mean( corrs );
}
void ExprPredictor::predict2(vector< vector< double > >& occs )
{
    ExprFunc* func = createExprFunc( par_model );
    vector< double > fOcc;
    vector< double > corrs;
    vector< vector< double > > occstemp(  seqSitesb.size(),vector< double >(seqSitesb[0].size()) ); // [0] is irrelevant since all experiments have same # of seqs
    for ( int i = 0; i < seqSitesb.size(); i++ ) {   // number of experiments
    for ( int j = 0; j < seqSitesb[ i ].size(); j++ ) {		// number of sequences
    vector< double > concs = factorExprData.getCol( 0 );  //103
    func->predictOcc( seqSitesb[ i ][ j ], i, concs,fOcc );
    occstemp[i][j] =  fOcc[i] ;
}
}
occs = occstemp;
}
double ExprPredictor::compAvgCorr2(  ExprPar& par )
{
    ExprFunc* func = createExprFunc( par );
    vector< double > corrs;
    Matrix f = factorExprData;
    Matrix e = exprData;
    double totalSim = 0;
    for ( int i = 0; i < 3; i++ ) {  // should be i < nSeqs()
    vector< double > predictedExprs;
    vector< double > observedExprs;
    for ( int j = 0; j < nConds(); j++ ) {
        vector< double > concs = factorExprData.getCol( j );
        if(j == 0 ) {
            anny.annoty2( seqsy[ i], seqSites[ i ], f, e, *func , j, ExprPredictor::seqNmes[i]);
        }
        else {anny.annoty2( seqsy[ i], seqSites[ i ], f, e, *func , j, " "); }
        double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
        predictedExprs.push_back( predicted );
        double observed = exprData( i, j );
        observedExprs.push_back( observed );
    }
    corrs.push_back( correlation( predictedExprs, observedExprs ) );
}
par_model = par;
cout <<  endl;
if( corrs != corrs) {
    cout << "corrs != corrs" << endl;
    return .5;
}
return mean( corrs );
}
double ExprPredictor::compAvgCorrborder8(  ExprPar& par )
{
    double expmin = .15;
    double expmax = .65;
    double bottomofdborder = .15;
    double topofdborder = .65;
    vector< int > initint;
    for ( int ab = 0; ab < nSeqs() ; ab++) {
        AllBorders.push_back( initint );
    }
    par_model = par;
    AllBorders.clear();
    AllData.clear();
    ExprFunc* func = createExprFunc( par );
    vector< double > corrs;
    vector< string > rowLabels(0);
    vector< string > colLabels(0);
    vector< vector< double > > data(0) ;
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * TT = gsl_rng_default;	// create rng type
    rng = gsl_rng_alloc( TT );
    gsl_rng_set( rng, 12345 );		// set the seed equal to simulTime(0)
    double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) {
        colLabels.push_back( label );
    }
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){
            vals.push_back( val  );
        }// while ss >> val
        data.push_back( vals );
    }
    Matrix data1( data );
    Matrix f = factorExprData;
    exprData = data;
    Matrix e = data;
    int nrow = exprData.nRows();
    int ncol =exprData.nCols();
    cout << " exprData nrows " << nrow << endl;
    gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    int tit;
    int tib;
    for (int i=0; i<nrow;i++) {
        vector< double > reD;
        reD = exprData.getRow(i);
        int mi;
        gsl_vector *rowexprData = gsl_vector_alloc(ncol);
        rowexprData = vector2gsl(reD);
        mi = gsl_vector_max_index(rowexprData);
        double m;
        m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
        double exptrace;
        double exppeek;
        tit = 0;
        tib =0;
        if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
        gsl_vector_set(transition_Indicest, i , tit);  // here ti = NULL
        gsl_vector_set(transition_Indicesb, i , tib);
        break;
    }
    else{
        for (int j=mi; j< ncol; j++)  {
            exptrace = gsl_vector_get(rowexprData,j);
            exppeek = gsl_vector_get(rowexprData,j+1);
            if ( exppeek > topofdborder ) {
                continue;
            }
            else {
                if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
                if(exppeek < topofdborder ) {
                    tit =j;
                    gsl_vector_set(transition_Indicest, i , tit);
                    int counter=0;
                    for(;;) {
                        exptrace = gsl_vector_get(rowexprData,tit+counter);
                        counter++;
                        exppeek = gsl_vector_get(rowexprData,tit + counter);
                        if ( exppeek > bottomofdborder ) {
                            continue;
                        }
                        else {
                            tib=tit + counter;
                            gsl_vector_set(transition_Indicesb, i , tib);
                            break;
                        }//else
                    } // for(;;)
                    int minindex;
                    minindex = gsl_vector_min_index(rowexprData);
                    if (tit == 0) {
                        tit = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                        tib = minindex;
                        gsl_vector_set(transition_Indicest, i , tit);
                        gsl_vector_set(transition_Indicesb, i , tib);
                    }
                    break;
                } // if < topofdborder
            }//  else
        }// for j
    }// else
}// for i
vector< double > predictedExprs;
vector< double > observedExprs;
double totalSim = 0;
double rms = 0;
cout << "transitionindicest "<< endl;
for (int i=0; i<nrow;i++) {
    cout << endl << gsl_vector_get(transition_Indicest,i) << endl;
}
cout << "transitionindicesb " << endl;
for (int i=0; i<nrow;i++) {
    cout << endl << gsl_vector_get(transition_Indicesb,i) << endl;
}
int c = 0;
for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
c++;
cell.push_back( gsl_vector_get(transition_Indicest,i) );
cell.push_back( gsl_vector_get(transition_Indicesb,i) );
if( gsl_vector_get(transition_Indicest,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    double observed = .7;  // .7 was for debugging, needs to change
    AllBorders[i].push_back( 0 );
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    AllBorders[i].push_back( 0 );
    double observed = .7;  // .7 was for debugging, needs to change
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
else{
    vector< double > concst = factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );
    anny.annotydorsal( seqsy[ i], seqSites[ i ], f, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
    AllBorders[i].push_back( gsl_vector_get(transition_Indicest,i) );
    double predictedt = func->predictExpr( seqSites[ i ], seqLengths[i], concst );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
    vector< double > concs2 = factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );
    anny.annotydorsal( seqsy[ i], seqSitesbot[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
    AllBorders[i].push_back( gsl_vector_get(transition_Indicesb,i) );
    double predicted2 = func->predictExpr( seqSitesbot[ i ], seqLengths[i], concs2 );
    predictedExprs.push_back( predicted2 );
    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
    observedExprs.push_back( observed2 );
    rms += ( predicted2 - observed2 )*( predicted2 - observed2 );   // this needs to be our target expression (on or off ?  or .5, shouldn't it be .5)
}
}  // for i
// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
bottomofdborder = .5;
topofdborder = .8;
nrow = factorExprData.nRows();
ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
cout << "transitionindicests "<< endl << gsl_vector_get(transition_Indicests,2) << endl;
cout << "transitionindicesbs " << endl <<gsl_vector_get(transition_Indicesbs,2) << endl;
for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.
vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );  // CHECK GENE i for nee, if not don't run annoty3.
if( concheck[2]!=0){ seqSitesm1[ i ]=seqSites[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
predictedExprs.push_back( predictedt );
double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
observedExprs.push_back( observedt );
rms += ( predictedt - observedt )*( predictedt - observedt );
continue; } // don't add more snail sites, if i is not an nee.
anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
predictedExprs.push_back( predictedt );
AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
observedExprs.push_back( observedt );
rms += ( predictedt - observedt )*( predictedt - observedt );
}
for ( int i = 0; i < nSeqs(); i++ ) {
    vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicesbs,2) );
    vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );  // CHECK GENE i for nee, if not don't run annoty3.
    if( concheck[2]!=0){seqSitesm2[ i ]=seqSitesbot[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
    double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
    continue; } // don't add more snail sites, if i is not an nee.
    anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
    double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
}
Matrix f2 =factorExprData;
vector< double > temp = factorExprData.getRow(0);
double noise = .1;
for( int i = 0 ; i < temp.size(); i++ ){
    temp[i]= temp[i] + noise;
    if(temp[i] > 1 ) { temp[i] = 1 ;}
}
f2.setRow(0,temp);
for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
if( gsl_vector_get(transition_Indicest,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    AllBorders[i].push_back( 0 );
    double observed = .7;  // .7 was for debugging, needs to change
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    AllBorders[i].push_back( 0 );
    double observed = .7;  // .7 was for debugging, needs to change
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
else{
    vector< double > concst = f2.getCol( gsl_vector_get(transition_Indicest,i) );
    anny.annotydorsal( seqsy[ i], seqSitesf2[ i ], f2, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
    AllBorders[i].push_back( gsl_vector_get(transition_Indicest,i) );
    double predictedt = func->predictExpr( seqSitesf2[ i ], seqLengths[i], concst );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
    vector< double > concs2 = f2.getCol( gsl_vector_get(transition_Indicesb,i) );
    anny.annotydorsal( seqsy[ i], seqSitesbotf2[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
    AllBorders[i].push_back( gsl_vector_get(transition_Indicesb,i) );
    double predicted2 = func->predictExpr( seqSitesbotf2[ i ], seqLengths[i], concs2 );
    predictedExprs.push_back( predicted2 );
    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
    observedExprs.push_back( observed2 );
    rms += ( predicted2 - observed2 )*( predicted2 - observed2 );   // this needs to be our target expression (on or off ?  or .5, shouldn't it be .5)
}
}  // for i
// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
bottomofdborder = .5;
topofdborder = .8;
nrow = factorExprData.nRows();
ncol =factorExprData.nCols();
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.
vector< double > concsm = f2.getCol( gsl_vector_get(transition_Indicests,2) );  // here we retrieve the snail border..
vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );  // CHECK GENE i for nee, if not don't run annoty3.
if( concheck[2]!=0){
    seqSitesm1f2[ i ]=seqSitesf2[i]; AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
    double predictedt = func->predictExpr( seqSitesm1f2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
    continue; } // don't add more snail sites, if i is not an nee.
    anny.annoty3( seqsy[ i], seqSitesm1f2[ i ], f2, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSitesf2[i]);
    double predictedt = func->predictExpr( seqSitesm1f2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
    double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
}
for ( int i = 0; i < nSeqs(); i++ ) {
    vector< double > concsm = f2.getCol( gsl_vector_get(transition_Indicesbs,2) );
    vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );  // CHECK GENE i for nee, if not don't run annoty3.
    if( concheck[2]!=0){ seqSitesm2f2[ i ]=seqSitesbotf2[i]; AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
    double predictedt = func->predictExpr( seqSitesm2f2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
    continue; } // don't add more snail sites, if i is not an nee.
    anny.annoty3( seqsy[ i], seqSitesm2f2[ i ], f2, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSitesf2[i]);
    AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
    double predictedt = func->predictExpr( seqSitesm2f2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
}
AllData.push_back(seqSites ) ;
AllData.push_back(seqSitesbot ) ;
AllData.push_back( seqSitesm1) ;
AllData.push_back(seqSitesm2 ) ;
AllData.push_back(seqSitesf2 ) ;
AllData.push_back(seqSitesbotf2 ) ;
AllData.push_back(seqSitesm1f2 ) ;
AllData.push_back(seqSitesm2f2 ) ;
obj_model = sqrt(rms/(nSeqs()*8));  // this was set for testing of correlation matrix of parameters, 10.5.11
return  obj_model;  // rms;
}
double ExprPredictor::compAvgCorrborder(  ExprPar& par )
{
    double expmin = .35;
    double expmax = .65;
    double bottomofdborder = .35;
    double topofdborder = .65;
    vector< int > initint;
    for ( int ab = 0; ab < nSeqs() ; ab++) {
        AllBorders.push_back( initint );
    }
    par_model = par;
    AllBorders.clear();
    AllData.clear();
    ExprFunc* func = createExprFunc( par );
    vector< double > corrs;
    vector< string > rowLabels(0);
    vector< string > colLabels(0);
    vector< vector< double > > data(0) ;
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * TT = gsl_rng_default;	// create rng type
    rng = gsl_rng_alloc( TT );
    gsl_rng_set( rng, 12345 );		// set the seed equal to simulTime(0)
    double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) {
        colLabels.push_back( label );
    }
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){
            vals.push_back( val  );
        }// while ss >> val
        data.push_back( vals );
    }
    Matrix data1( data );
    Matrix f = factorExprData;
    exprData = data;
    Matrix e = data;
    int nrow = exprData.nRows();
    int ncol =exprData.nCols();
    gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    int tit;
    int tib;
    for (int i=0; i<nrow;i++) {
        vector< double > reD;
        reD = exprData.getRow(i);
        int mi;
        gsl_vector *rowexprData = gsl_vector_alloc(ncol);
        rowexprData = vector2gsl(reD);
        mi = gsl_vector_max_index(rowexprData);
        double m;
        m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
        double exptrace;
        double exppeek;
        tit = 0;
        tib =0;
        if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
        gsl_vector_set(transition_Indicest, i , tit);  // here ti = NULL
        gsl_vector_set(transition_Indicesb, i , tib);
        break;
    }
    else{
        for (int j=mi; j< ncol; j++)  {
            exptrace = gsl_vector_get(rowexprData,j);
            exppeek = gsl_vector_get(rowexprData,j+1);
            if ( exppeek > topofdborder ) {
                continue;
            }
            else {
                if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
                if(exppeek < topofdborder ) {
                    tit =j;
                    gsl_vector_set(transition_Indicest, i , tit);
                    int counter=0;
                    for(;;) {
                        exptrace = gsl_vector_get(rowexprData,tit+counter);
                        counter++;
                        exppeek = gsl_vector_get(rowexprData,tit + counter);
                        if ( exppeek > bottomofdborder ) {
                            continue;
                        }
                        else {
                            tib=tit + counter;
                            gsl_vector_set(transition_Indicesb, i , tib);
                            break;
                        }//else
                    } // for(;;)
                    int minindex;
                    minindex = gsl_vector_min_index(rowexprData);
                    if (tit == 0) {
                        tit = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                        tib = minindex;
                        gsl_vector_set(transition_Indicest, i , tit);
                        gsl_vector_set(transition_Indicesb, i , tib);
                    }
                    break;
                } // if < topofdborder
            }//  else
        }// for j
    }// else
}// for i
vector< double > predictedExprs;
vector< double > observedExprs;
double totalSim = 0;
double rms = 0;
for (int i=0; i<nrow;i++) {
}
for (int i=0; i<nrow;i++) {
}
int c = 0;
for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
c++;
cell.push_back( gsl_vector_get(transition_Indicest,i) );
cell.push_back( gsl_vector_get(transition_Indicesb,i) );
if( gsl_vector_get(transition_Indicest,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    double observed = .7;  // .7 was for debugging, needs to change
    AllBorders[i].push_back( 0 );
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    AllBorders[i].push_back( 0 );
    double observed = .7;  // .7 was for debugging, needs to change
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
else{
    vector< double > concst = factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );
    anny.annotydorsal( seqsy[ i], seqSites[ i ], f, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
    AllBorders[i].push_back( gsl_vector_get(transition_Indicest,i) );
    double predictedt = func->predictExpr( seqSites[ i ], seqLengths[i], concst );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
    vector< double > concs2 = factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );
    anny.annotydorsal( seqsy[ i], seqSitesbot[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
    AllBorders[i].push_back( gsl_vector_get(transition_Indicesb,i) );
    double predicted2 = func->predictExpr( seqSitesbot[ i ], seqLengths[i], concs2 );
    predictedExprs.push_back( predicted2 );
    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
    observedExprs.push_back( observed2 );
    rms += ( predicted2 - observed2 )*( predicted2 - observed2 );   // this needs to be our target expression (on or off ?  or .5, shouldn't it be .5)
}
}  // for i
// remember the transition indices start at 0, while the spreadsheet labels the columns starting at 1.
bottomofdborder = .5;
topofdborder = .8;
nrow = factorExprData.nRows();
ncol =factorExprData.nCols();
gsl_vector *transition_Indicesbs = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
gsl_vector *transition_Indicests = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
int tits;
int tibs;
for (int i=0; i<nrow;i++) {
    vector< double > reD;
    reD = factorExprData.getRow(i);
    int mi;
    gsl_vector *rowexprData = gsl_vector_alloc(ncol);
    rowexprData = vector2gsl(reD);
    mi = gsl_vector_max_index(rowexprData);
    double m;
    m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
    double exptrace;
    double exppeek;
    tits = 0;
    tibs =0;
    if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
    gsl_vector_set(transition_Indicests, i , tits);  // here ti = NULL
    gsl_vector_set(transition_Indicesbs, i , tibs);
    break;
}
else{
    for (int j=mi; j< ncol; j++)  {
        exptrace = gsl_vector_get(rowexprData,j);
        exppeek = gsl_vector_get(rowexprData,j+1);
        if ( exppeek > topofdborder ) {
            continue;
        }
        else {
            if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
            if(exppeek < topofdborder ) {
                tits =j;
                gsl_vector_set(transition_Indicests, i , tits);
                int counter=0;
                for(;;) {
                    exptrace = gsl_vector_get(rowexprData,tits+counter);
                    counter++;
                    exppeek = gsl_vector_get(rowexprData,tits + counter);
                    if ( exppeek > bottomofdborder ) {
                        continue;
                    }
                    else {
                        tibs=tits + counter;
                        gsl_vector_set(transition_Indicesbs, i , tibs);
                        break;
                    }//else
                } // for(;;)
                int minindex;
                minindex = gsl_vector_min_index(rowexprData);
                if (tits == 0) {
                    tits = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                    tibs = minindex;
                    gsl_vector_set(transition_Indicests, i , tits);
                    gsl_vector_set(transition_Indicesbs, i , tibs);
                }
                break;
            } // if < topofdborder
        }//  else
    }// for j
}// else
}//for i
for ( int i = 0; i < nSeqs(); i++ ) {  // given the neuroectoderm structures, we calculate (estimate?) the mesoderms structure.
vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicests,2) );
vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicest,i) );  // CHECK GENE i for nee, if not don't run annoty3.
if( concheck[2]!=0){ seqSitesm1[ i ]=seqSites[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
predictedExprs.push_back( predictedt );
double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
observedExprs.push_back( observedt );
rms += ( predictedt - observedt )*( predictedt - observedt );
cout << " rms " << rms << endl;
continue; } // don't add more snail sites, if i is not an nee.
anny.annoty3( seqsy[ i], seqSitesm1[ i ], f, e, *func , gsl_vector_get(transition_Indicests,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
double predictedt = func->predictExpr( seqSitesm1[ i ], seqLengths[i], concsm );
predictedExprs.push_back( predictedt );
AllBorders[i].push_back( gsl_vector_get(transition_Indicests,2) );
double observedt = exprData( i, gsl_vector_get(transition_Indicests,2) );
observedExprs.push_back( observedt );
rms += ( predictedt - observedt )*( predictedt - observedt );
}
for ( int i = 0; i < nSeqs(); i++ ) {
    vector< double > concsm = factorExprData.getCol( gsl_vector_get(transition_Indicesbs,2) );
    vector<double > concheck=factorExprData.getCol( gsl_vector_get(transition_Indicesb,i) );  // CHECK GENE i for nee, if not don't run annoty3.
    if( concheck[2]!=0){seqSitesm2[ i ]=seqSitesbot[i];AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
    double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
    continue; } // don't add more snail sites, if i is not an nee.
    anny.annoty3( seqsy[ i], seqSitesm2[ i ], f, e, *func , gsl_vector_get(transition_Indicesbs,2), ExprPredictor::seqNmes[i],i, seqSites[i]);
    double predictedt = func->predictExpr( seqSitesm2[ i ], seqLengths[i], concsm );
    predictedExprs.push_back( predictedt );
    AllBorders[i].push_back( gsl_vector_get(transition_Indicesbs,2) );
    double observedt = exprData( i, gsl_vector_get(transition_Indicesbs,2));
    observedExprs.push_back( observedt );
    rms += ( predictedt - observedt )*( predictedt - observedt );
}
AllData.push_back(seqSites ) ;
AllData.push_back(seqSitesbot ) ;
AllData.push_back( seqSitesm1) ;
AllData.push_back(seqSitesm2 ) ;
obj_model = sqrt(rms/(nSeqs()*4));  // this was set for testing of correlation matrix of parameters, 10.5.11
return  obj_model;  // rms;
}
double ExprPredictor::compAvgCorrborder2(  ExprPar& par )
{
    double expmin = .35;
    double expmax = .65;
    double bottomofdborder = .35;
    double topofdborder = .65;
    par_model = par;
    ExprFunc* func = createExprFunc( par );
    vector< double > corrs;
    vector< string > rowLabels(0);
    vector< string > colLabels(0);
    vector< vector< double > > data(0) ;
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * TT = gsl_rng_default;	// create rng type
    rng = gsl_rng_alloc( TT );
    gsl_rng_set( rng, 12345 );		// set the seed equal to simulTime(0)
    double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) {
        colLabels.push_back( label );
    }
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){
            vals.push_back( val  );
        }// while ss >> val
        data.push_back( vals );
    }
    Matrix data1( data );
    Matrix f = factorExprData;
    exprData = data;
    Matrix e = data;
    int nrow = exprData.nRows();
    int ncol =exprData.nCols();
    gsl_vector *transition_Indicesb = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    gsl_vector *transition_Indicest = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    int tit;
    int tib;
    for (int i=0; i<nrow;i++) {
        vector< double > reD;
        reD = exprData.getRow(i);
        int mi;
        gsl_vector *rowexprData = gsl_vector_alloc(ncol);
        rowexprData = vector2gsl(reD);
        mi = gsl_vector_max_index(rowexprData);
        double m;
        m = gsl_vector_max(rowexprData);      // the max function starts from the left and works its way to the right (starts with lowest indices)
        double exptrace;
        double exppeek;
        tit = 0;
        tib =0;
        if( m < bottomofdborder ) {                           // this is a security check to make sure borders exist
        gsl_vector_set(transition_Indicest, i , tit);  // here ti = NULL
        gsl_vector_set(transition_Indicesb, i , tib);
        break;
    }
    else{
        for (int j=mi; j< ncol; j++)  {
            exptrace = gsl_vector_get(rowexprData,j);
            exppeek = gsl_vector_get(rowexprData,j+1);
            if ( exppeek > topofdborder ) {
                continue;
            }
            else {
                if( exptrace == exppeek) { continue; }  // peek ahead to make sure trace is not on a saddle point (plateu)
                if(exppeek < topofdborder ) {
                    tit =j;
                    gsl_vector_set(transition_Indicest, i , tit);
                    int counter=0;
                    for(;;) {
                        exptrace = gsl_vector_get(rowexprData,tit+counter);
                        counter++;
                        exppeek = gsl_vector_get(rowexprData,tit + counter);
                        if ( exppeek > bottomofdborder ) {
                            continue;
                        }
                        else {
                            tib=tit + counter;
                            gsl_vector_set(transition_Indicesb, i , tib);
                            break;
                        }//else
                    } // for(;;)
                    int minindex;
                    minindex = gsl_vector_min_index(rowexprData);
                    if (tit == 0) {
                        tit = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
                        tib = minindex;
                        gsl_vector_set(transition_Indicest, i , tit);
                        gsl_vector_set(transition_Indicesb, i , tib);
                    }
                    break;
                } // if < topofdborder
            }//  else
        }// for j
    }// else
}// for i
vector< double > predictedExprs;
vector< double > observedExprs;
double totalSim = 0;
double rms = 0;
Matrix f2 =factorExprData;
vector< double > temp = factorExprData.getRow(0);
double noise = .1;
for( int i = 0 ; i < temp.size(); i++ ){
    temp[i]= temp[i] + noise;
    if(temp[i] > 1 ) { temp[i] = 1 ;}
}
for ( int i = 0; i < nSeqs(); i++ ) {  // should be i < nSeqs()
if( gsl_vector_get(transition_Indicest,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    cout << " gsl_transindex = 0 in f2 " << endl;
    double observed = .7;  // .7 was for debugging, needs to change
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
if( gsl_vector_get(transition_Indicesb,i) == 0 ) {
    double predicted = 0;   // func->predictExpr( seqSites[ i ], seqLengths[i], concs );
    predictedExprs.push_back( predicted );
    cout << " gslindex = 0 " << endl;
    double observed = .7;  // .7 was for debugging, needs to change
    observedExprs.push_back( observed );
    rms += abs( predicted - 0 );
}
else{
    vector< double > concst = f2.getCol( gsl_vector_get(transition_Indicest,i) );
    anny.annotydorsal( seqsy[ i], seqSitesf2[ i ], f2, e, *func , gsl_vector_get(transition_Indicest,i), ExprPredictor::seqNmes[i]);
    double predictedt = func->predictExpr( seqSitesf2[ i ], seqLengths[i], concst );
    predictedExprs.push_back( predictedt );
    double observedt = exprData( i, gsl_vector_get(transition_Indicest,i) );
    observedExprs.push_back( observedt );
    rms += abs( predictedt - observedt );
    vector< double > concs2 = f2.getCol( gsl_vector_get(transition_Indicesb,i) );
    anny.annotydorsal( seqsy[ i], seqSitesbotf2[ i ], f, e, *func , gsl_vector_get(transition_Indicesb,i), ExprPredictor::seqNmes[i]);
    double predicted2 = func->predictExpr( seqSitesbotf2[ i ], seqLengths[i], concs2 );
    predictedExprs.push_back( predicted2 );
    double observed2 = exprData( i, gsl_vector_get(transition_Indicesb,i) );
    observedExprs.push_back( observed2 );
}
}  // for i
return rms/nSeqs();
}
int Random(int n)
{
    return rand() % n ;
}
int ExprPredictor::createccdata()
{
    double expmin = .35;
    double expmax = .65;
    int k=3;
    vector< string > rowLabels(0);
    vector< string > colLabels(0);
    vector< vector< double > > data(0) ;
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * TT = gsl_rng_default;	// create rng type
    rng = gsl_rng_alloc( TT );
    gsl_rng_set( rng, 12345 );		// set the seed equal to simulTime(0)
    double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) colLabels.push_back( label );
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){
            rand_site_index = gsl_ran_flat( rng, -.15, .15 );
            if( val + rand_site_index > 1 ) {
                vals.push_back( 1 );
            }
            if( val + rand_site_index < 0 ) {
                vals.push_back( 0 );
            }
            if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
                vals.push_back( val + rand_site_index );
            }
        } // while ss >> val
        data.push_back( vals );
    }
    Matrix data1( data );
    Matrix f = factorExprData;
    exprData = data;
    Matrix e = data;
    int nrow = exprData.nRows();
    int ncol =exprData.nCols();
    gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    for (int i=0; i<nrow;i++) {
        vector< double > reD;
        reD = exprData.getRow(i);
        int mi;
        gsl_vector *rowexprData = gsl_vector_alloc(ncol);
        rowexprData = vector2gsl(reD);
        mi = gsl_vector_max_index(rowexprData);
        double m;
        m = gsl_vector_max(rowexprData);
        int ti;
        ti = 0;
        if( m < expmin ) {
            gsl_vector_set(transition_Indices, i , ti);  // here ti = NULL
            break;
        }
        else{
            for (int j=mi; j< ncol; j++)  {
                if (expmax < gsl_vector_get(rowexprData,j) ) {
                    continue;
                }
                else {
                    ti = j; break;
                }
            }
            int minindex;
            minindex = gsl_vector_min_index(rowexprData);
            if (ti == 0) {
                ti = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
            }
            gsl_vector_set(transition_Indices, i , ti);                    // pushback the ti that is the transition index..
        }// else
    }// for i
    for (int i=0; i<nrow;i++) {
        cout << "trans index : " << gsl_vector_get(transition_Indices,i) << endl;
    }
    ifstream seq_file("efastanotw10txt");  //("efastanotw6.txt");
    ifstream expr_file( "expre10.tab");//"expre6.tab");
    assert (seq_file.is_open() &&  expr_file.is_open());
    string header;
    getline(expr_file, header);
    string temp;
    vector <string> seq_name;
    vector <string> seq;
    vector <string> expr;
    seq.clear();
    expr.clear();
    while(!seq_file.eof()){
        temp = "";
        getline(seq_file, temp);
        if(temp.length() == 0){
            break;
        }
        string name1 (temp, 1, temp.length() - 1);
        seq_name.push_back(name1);
        getline(seq_file, temp);
        seq.push_back(temp);
        getline(expr_file, temp);
        expr.push_back(temp);
    }
    int size = seq.size();
    ofstream meso_seq_file ("mesoseq.txt");
    ofstream meso_expr_file ("mesoexpr.txt");
    ofstream neuro_seq_file("neuroseq.txt");
    ofstream neuro_expr_file("neuroexpr.txt");
    assert (meso_seq_file.is_open() &&  meso_expr_file.is_open() && neuro_seq_file.is_open() && neuro_expr_file.is_open());
    meso_expr_file << header << endl;
    neuro_expr_file << header << endl;
    for(int i = 0; i < size; i++){
        if( gsl_vector_get(transition_Indices,i) < 10 ) {
            meso_seq_file << ">" << seq_name[i] << endl;
            meso_seq_file << seq[i] << endl;
            meso_expr_file << expr[i] << endl;
        }
        else{
            neuro_seq_file << ">" << seq_name[i] << endl;
            neuro_seq_file << seq[i] << endl;
            neuro_expr_file << expr[i] << endl;
        }
    }
    neuro_seq_file.close();
    neuro_expr_file.close();
    meso_seq_file.close();
    meso_expr_file.close();
    seq_file.close();
    expr_file.close();
    ifstream seq_filem("mesoseq.txt");
    ifstream expr_filem("mesoexpr.txt");
    string headerm;
    getline(expr_filem, headerm);
    string tempm;
    vector <string> seq_namem;
    vector <string> seqm;
    vector <string> exprm;
    seqm.clear();
    exprm.clear();
    while(!seq_filem.eof()){
        tempm = "";
        getline(seq_filem, tempm);
        if(tempm.length() == 0){
            break;
        }
        string namem (tempm, 1, tempm.length() - 1);
        seq_namem.push_back(namem);
        getline(seq_filem, tempm);
        seqm.push_back(tempm);
        getline(expr_filem, tempm);
        exprm.push_back(tempm);
    }
    int sizem = seqm.size();
    for(int i = 0; i < seqm.size(); i++){
    }
    vector <int> indicesm(sizem, 0);
    for(int i = 0; i < indicesm.size(); i++)
    indicesm[i] = i;
    std::srand(12345);
    random_shuffle(indicesm.begin(), indicesm.end(), Random);
    int min_sizem = sizem/k;
    vector <int> partition_sizesm(k, min_sizem);
    int residuem = sizem - k * min_sizem;
    for(int i = 0; i < residuem; i++){
        partition_sizesm[i] ++;
    }
    int summ = 0;
    for(int i = 0; i < k; i++){
        summ += partition_sizesm[i];
    }
    assert (summ == sizem);
    int startm = 0;
    vector <int> startsm;
    startsm.clear();
    for(int i = 0; i < k; i++){
        startsm.push_back(startm);
        startm += partition_sizesm[i];
    }
    for(int i = 0; i < k; i++){
    }
    for(int i = 0; i < k; i++){
        ostringstream numberm;
        numberm << i;
        ofstream train_seq_file (("train_seq_" + numberm.str()).c_str());
        ofstream train_expr_file (("train_expr_" + numberm.str()).c_str());
        ofstream test_seq_file(("test_seq_" + numberm.str()).c_str());
        ofstream test_expr_file(("test_expr_" + numberm.str()).c_str());
        assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
        train_expr_file << headerm << endl;
        test_expr_file << headerm << endl;
        for(int j = 0; j < sizem; j++){
            int indexm = indicesm[j];
            if(j >= startsm[i] && j < startsm[i] + partition_sizesm[i]){
                test_seq_file << ">" << seq_namem[indexm] << endl;
                test_seq_file << seqm[indexm] << endl;
                test_expr_file << exprm[indexm] << endl;
            }
            else{
                train_seq_file << ">" << seq_namem[indexm] << endl;
                train_seq_file << seqm[indexm] << endl;
                train_expr_file << exprm[indexm] << endl;
            }
        }
        train_seq_file.close();
        train_expr_file.close();
        test_seq_file.close();
        test_expr_file.close();
    }
    ifstream seq_filen("neuroseq.txt");
    ifstream expr_filen("neuroexpr.txt");
    string headern;
    getline(expr_filen, headern);
    string tempn;
    vector <string> seq_namen;
    vector <string> seqn;
    vector <string> exprn;
    seqn.clear();
    exprn.clear();
    while(!seq_filen.eof()){
        tempn = "";
        getline(seq_filen, tempn);
        if(tempn.length() == 0){
            break;
        }
        string namen (tempn, 1, tempn.length() - 1);
        seq_namen.push_back(namen);
        getline(seq_filen, tempn);
        seqn.push_back(tempn);
        getline(expr_filen, tempn);
        exprn.push_back(tempn);
    }
    int sizen = seqn.size();
    for(int i = 0; i < seqn.size(); i++){
    }
    vector <int> indicesn(sizen, 0);
    for(int i = 0; i < indicesn.size(); i++)
    indicesn[i] = i;
    std::srand(12345);
    random_shuffle(indicesn.begin(), indicesn.end(), Random);
    int min_sizen = sizen/k;
    vector <int> partition_sizesn(k, min_sizen);
    int residuen = sizen - k * min_sizen;
    for(int i = 0; i < residuen; i++){
        partition_sizesn[i] ++;
    }
    int sumn = 0;
    for(int i = 0; i < k; i++){
        sumn += partition_sizesn[i];
    }
    assert (sumn == sizen);
    int startn = 0;
    vector <int> startsn;
    startsn.clear();
    for(int i = 0; i < k; i++){
        startsn.push_back(startn);
        startn += partition_sizesn[i];
    }
    for(int i = 0; i < k; i++){
    }
    for(int i = 0; i < k; i++){
        ostringstream numbern;
        numbern << i;
        ofstream train_seq_file (("train_seq_" + numbern.str()).c_str(), ios::app);
        ofstream train_expr_file (("train_expr_" + numbern.str()).c_str(), ios::app);
        ofstream test_seq_file(("test_seq_" + numbern.str()).c_str(), ios::app);
        ofstream test_expr_file(("test_expr_" + numbern.str()).c_str(), ios::app);
        assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
        for(int j = 0; j < sizen; j++){
            int indexn = indicesn[j];
            if(j >= startsn[i] && j < startsn[i] + partition_sizesn[i]){
                test_seq_file << ">" << seq_namen[indexn] << endl;
                test_seq_file << seqn[indexn] << endl;
                test_expr_file << exprn[indexn] << endl;
            }
            else{
                train_seq_file<< ">" << seq_namen[indexn] << endl;
                train_seq_file << seqn[indexn] << endl;
                train_expr_file << exprn[indexn] << endl;
            }
        }
        train_seq_file.close();
        train_expr_file.close();
        test_seq_file.close();
        test_expr_file.close();
    }
    return 1;
}
int ExprPredictor::createccdata2()
{
    double expmin = .35;
    double expmax = .65;
    int k=3;
    vector< string > rowLabels(0);
    vector< string > colLabels(0);
    vector< vector< double > > data(0) ;
    ifstream fin( file.c_str() );
    if ( !fin ) {
        return RET_ERROR;
    }
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * TT = gsl_rng_default;	// create rng type
    rng = gsl_rng_alloc( TT );
    gsl_rng_set( rng, 12345 );		// set the seed equal to simulTime(0)
    double  rand_site_index;	// from ~/C++exercises/Bins/Mscan/wtmx_scanmc1116.cpp
    rowLabels.clear();
    colLabels.clear();
    data.clear();
    string line, first, label;
    getline( fin, line );
    stringstream ss( line );
    ss >> first;
    while ( ss >> label ) colLabels.push_back( label );
    while ( !fin.eof() ) {
        string line;
        getline( fin, line );
        if ( line.empty() ) continue;
        stringstream ss( line );
        string name;
        ss >> name;
        rowLabels.push_back( name );
        double val;
        vector< double > vals;
        while ( ss >> val ){
            rand_site_index = gsl_ran_flat( rng, -.15, .15 );
            if( val + rand_site_index > 1 ) {
                vals.push_back( 1 );
            }
            if( val + rand_site_index < 0 ) {
                vals.push_back( 0 );
            }
            if( val + rand_site_index >= 0  && val + rand_site_index <= 1 ) {
                vals.push_back( val + rand_site_index );
            }
        } // while ss >> val
        data.push_back( vals );
    }
    Matrix data1( data );
    Matrix f = factorExprData;
    exprData = data;
    Matrix e = data;
    int nrow = exprData.nRows();
    int ncol =exprData.nCols();
    gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    for (int i=0; i<nrow;i++) {
        vector< double > reD;
        reD = exprData.getRow(i);
        int mi;
        gsl_vector *rowexprData = gsl_vector_alloc(ncol);
        rowexprData = vector2gsl(reD);
        mi = gsl_vector_max_index(rowexprData);
        double m;
        m = gsl_vector_max(rowexprData);
        int ti;
        ti = 0;
        if( m < expmin ) {
            gsl_vector_set(transition_Indices, i , ti);  // here ti = NULL
            break;
        }
        else{
            for (int j=mi; j< ncol; j++)  {
                if (expmax < gsl_vector_get(rowexprData,j) ) {
                    continue;
                }
                else {
                    ti = j; break;
                }
            }
            int minindex;
            minindex = gsl_vector_min_index(rowexprData);
            if (ti == 0) {
                ti = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
            }
            gsl_vector_set(transition_Indices, i , ti);                    // pushback the ti that is the transition index..
        }// else
    }// for i
    for (int i=0; i<nrow;i++) {
    }
    ifstream seq_file("efastanotw10.txt");  //("efastanotw6.txt");
    ifstream expr_file( "expre10.tab");//"expre6.tab");
    assert (seq_file.is_open() &&  expr_file.is_open());
    string header;
    getline(expr_file, header);
    string temp;
    vector <string> seq_name;
    vector <string> seq;
    vector <string> expr;
    seq.clear();
    expr.clear();
    while(!seq_file.eof()){
        temp = "";
        getline(seq_file, temp);
        if(temp.length() == 0){
            break;
        }
        string name1 (temp, 1, temp.length() - 1);
        seq_name.push_back(name1);
        getline(seq_file, temp);
        seq.push_back(temp);
        getline(expr_file, temp);
        expr.push_back(temp);
    }
    int size = seq.size();
    ofstream meso_seq_file ("mesoseq.txt");
    ofstream meso_expr_file ("mesoexpr.txt");
    ofstream neuro_seq_file("neuroseq.txt");
    ofstream neuro_expr_file("neuroexpr.txt");
    assert (meso_seq_file.is_open() &&  meso_expr_file.is_open() && neuro_seq_file.is_open() && neuro_expr_file.is_open());
    meso_expr_file << header << endl;
    neuro_expr_file << header << endl;
    for(int i = 0; i < size; i++){
        ostringstream numberm;
        numberm << i;
        if( gsl_vector_get(transition_Indices,i) < 10 ) {
            meso_seq_file << ">" << seq_name[i] << endl;
            meso_seq_file << seq[i] << endl;
            meso_expr_file << expr[i] << endl;
        }
        else{
            neuro_seq_file << ">" << seq_name[i] << endl;
            neuro_seq_file << seq[i] << endl;
            neuro_expr_file << expr[i] << endl;
        }
    }
    neuro_seq_file.close();
    neuro_expr_file.close();
    meso_seq_file.close();
    meso_expr_file.close();
    seq_file.close();
    expr_file.close();
    ifstream seq_filem("mesoseq.txt");
    ifstream expr_filem("mesoexpr.txt");
    string headerm;
    getline(expr_filem, headerm);
    string tempm;
    vector <string> seq_namem;
    vector <string> seqm;
    vector <string> exprm;
    seqm.clear();
    exprm.clear();
    while(!seq_filem.eof()){
        tempm = "";
        getline(seq_filem, tempm);
        if(tempm.length() == 0){
            break;
        }
        string namem (tempm, 1, tempm.length() - 1);
        seq_namem.push_back(namem);
        getline(seq_filem, tempm);
        seqm.push_back(tempm);
        getline(expr_filem, tempm);
        exprm.push_back(tempm);
    }
    int sizem = seqm.size();
    for(int i = 0; i < seqm.size(); i++){
    }
    vector <int> indicesm(sizem, 0);
    for(int i = 0; i < indicesm.size(); i++)
    indicesm[i] = i;
    std::srand(12345);
    random_shuffle(indicesm.begin(), indicesm.end(), Random);
    int min_sizem = sizem/k;
    vector <int> partition_sizesm(k, min_sizem);
    int residuem = sizem - k * min_sizem;
    for(int i = 0; i < residuem; i++){
        partition_sizesm[i] ++;
    }
    int summ = 0;
    for(int i = 0; i < k; i++){
        summ += partition_sizesm[i];
    }
    assert (summ == sizem);
    int startm = 0;
    vector <int> startsm;
    startsm.clear();
    for(int i = 0; i < k; i++){
        startsm.push_back(startm);
        startm += partition_sizesm[i];
    }
    for(int i = 0; i < k; i++){
    }
    for(int i = 0; i < k; i++){
        ostringstream numberm;
        numberm << i;
        ofstream train_seq_file (("train_seq_" + numberm.str()).c_str());
        ofstream train_expr_file (("train_expr_" + numberm.str()).c_str());
        ofstream test_seq_file(("test_seq_" + numberm.str()).c_str());
        ofstream test_expr_file(("test_expr_" + numberm.str()).c_str());
        assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
        train_expr_file << headerm << endl;
        test_expr_file << headerm << endl;
        for(int j = 0; j < sizem; j++){
            int indexm = indicesm[j];
            if(j >= startsm[i] && j < startsm[i] + partition_sizesm[i]){
                test_seq_file << ">" << seq_namem[indexm] << endl;
                test_seq_file << seqm[indexm] << endl;
                test_expr_file << exprm[indexm] << endl;
            }
            else{
                train_seq_file << ">" << seq_namem[indexm] << endl;
                train_seq_file << seqm[indexm] << endl;
                train_expr_file << exprm[indexm] << endl;
            }
        }
        train_seq_file.close();
        train_expr_file.close();
        test_seq_file.close();
        test_expr_file.close();
    }
    ifstream seq_filen("neuroseq.txt");
    ifstream expr_filen("neuroexpr.txt");
    string headern;
    getline(expr_filen, headern);
    string tempn;
    vector <string> seq_namen;
    vector <string> seqn;
    vector <string> exprn;
    seqn.clear();
    exprn.clear();
    while(!seq_filen.eof()){
        tempn = "";
        getline(seq_filen, tempn);
        if(tempn.length() == 0){
            break;
        }
        string namen (tempn, 1, tempn.length() - 1);
        seq_namen.push_back(namen);
        getline(seq_filen, tempn);
        seqn.push_back(tempn);
        getline(expr_filen, tempn);
        exprn.push_back(tempn);
    }
    int sizen = seqn.size();
    for(int i = 0; i < seqn.size(); i++){
    }
    vector <int> indicesn(sizen, 0);
    for(int i = 0; i < indicesn.size(); i++)
    indicesn[i] = i;
    std::srand(12345);
    random_shuffle(indicesn.begin(), indicesn.end(), Random);
    int min_sizen = sizen/k;
    vector <int> partition_sizesn(k, min_sizen);
    int residuen = sizen - k * min_sizen;
    for(int i = 0; i < residuen; i++){
        partition_sizesn[i] ++;
    }
    int sumn = 0;
    for(int i = 0; i < k; i++){
        sumn += partition_sizesn[i];
    }
    assert (sumn == sizen);
    int startn = 0;
    vector <int> startsn;
    startsn.clear();
    for(int i = 0; i < k; i++){
        startsn.push_back(startn);
        startn += partition_sizesn[i];
    }
    for(int i = 0; i < k; i++){
    }
    for(int i = 0; i < k; i++){
        ostringstream numbern;
        numbern << i;
        ofstream train_seq_file (("train_seq_" + numbern.str()).c_str(), ios::app);
        ofstream train_expr_file (("train_expr_" + numbern.str()).c_str(), ios::app);
        ofstream test_seq_file(("test_seq_" + numbern.str()).c_str(), ios::app);
        ofstream test_expr_file(("test_expr_" + numbern.str()).c_str(), ios::app);
        assert (train_seq_file.is_open() && train_expr_file.is_open() && test_seq_file.is_open() && test_expr_file.is_open());
        for(int j = 0; j < sizen; j++){
            int indexn = indicesn[j];
            if(j >= startsn[i] && j < startsn[i] + partition_sizesn[i]){
                test_seq_file << ">" << seq_namen[indexn] << endl;
                test_seq_file << seqn[indexn] << endl;
                test_expr_file << exprn[indexn] << endl;
            }
            else{
                train_seq_file << ">" << seq_namen[indexn] << endl;
                train_seq_file << seqn[indexn] << endl;
                train_expr_file << exprn[indexn] << endl;
            }
        }
        train_seq_file.close();
        train_expr_file.close();
        test_seq_file.close();
        test_expr_file.close();
    }
    return 1;
}
void ExprPredictor::compvar(vector< double >& vars) // added from Binsrep folder check READMEd there
{
    vector< double > predictedExprs;
    vector< double > stepobservedExprs;
    ExprFunc* func = createExprFunc( par_model );
    double rss = 0;
    double step = .001;
    vector< double > concs = factorExprData.getCol( 0 );
    vector< double > concs2 = concs;
    concs2[0] = concs[0] - step;
    for ( int i = 0; i < nSeqs(); i++ ) {
        double predicted = func->predictExpr( seqSites[ i ], 5, concs );
        predictedExprs.push_back( predicted );
        double stepobserved = func->predictExpr(seqSites[ i ], 5, concs2 );       //exprData( i, 0 );
        stepobservedExprs.push_back( stepobserved );
    }
    double partialcu = .5 * ( concs2[0] + concs[0] );
    assert( predictedExprs.size() == stepobservedExprs.size());
    for ( int i = 0; i <  predictedExprs.size(); i++ ) {
        vars.push_back(  partialcu * (predictedExprs[i]  -  stepobservedExprs[i]  ) / step )   ;
    }
}
double ExprPredictor::compAvgCrossCorr( const ExprPar& par ) const
{
    ExprFunc* func = createExprFunc( par );
    double totalSim = 0;
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > predictedExprs;
        vector< double > observedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
            vector< double > concs = factorExprData.getCol( j );
            double predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs );
            predictedExprs.push_back( predicted );
            double observed = exprData( i, j );
            observedExprs.push_back( observed );
        }
        totalSim += exprSimCrossCorr( predictedExprs, observedExprs );
    }
    return totalSim / nSeqs();
}
void ExprPredictor::compOccMat(const gsl_vector* v, void* params  )
{
    int nrow = exprData.nRows();
    int ncol =exprData.nCols();
    gsl_vector *transition_Indices = gsl_vector_alloc(nrow);  // vector that hold the ti for each gene..
    for (int i=0; i<nrow;i++) {
        vector< double > reD;
        reD = exprData.getRow(i);
        int mi;
        gsl_vector *rowexprData = gsl_vector_alloc(ncol);
        rowexprData = vector2gsl(reD);
        mi = gsl_vector_max_index(rowexprData);
        double m;
        m = gsl_vector_max(rowexprData);
        int ti;
        ti = 0;
        for (int j=mi; j< ncol-mi; j++)  {
            double temp;
            temp = m/2;
            if (temp < gsl_vector_get(rowexprData,j) ) {
                continue;
            }
            else {
                ti = j; break;
            }
        }
        int minindex;
        minindex = gsl_vector_min_index(rowexprData);
        if (ti == 0) {
            ti = minindex;            //some sequences are all zero or all the same value, which causes segemtation fault
        }
        gsl_vector_set(transition_Indices, i , ti);                    // pushback the ti that is the transition index..
    }
    ExprFunc* funcmat = createExprFunc2( par_model );
    gsl_matrix * Occupancy = gsl_matrix_alloc(nrow,par_model.nFactors() );  //   this should be  done in the optimization do loop
    for ( int i = 0; i < nSeqs(); i++ ) {
        vector< double > concs = factorExprData.getCol( gsl_vector_get(transition_Indices,i) );
        gsl_vector * Occvector =gsl_vector_alloc(motifs.size()); // there are 4 w's
        vector <double > fOcc(3);
        funcmat->predictOcc( seqSites[ i ], seqLengths[i], concs, fOcc); // don't need TranI or i, since concs only has the transition
        Occvector = vector2gsl(fOcc);  // for the fix_pars that don't have  an occupancy we need to make sure factorOcc has 1 initialized
        gsl_matrix_set_row(Occupancy , i , Occvector);
    }
    int fixedsize = fix_pars.size();
    int rows = Occupancy -> size1;
    int columns = Occupancy -> size2;
    int i,j,k;
    k=0;
    gsl_matrix *X = gsl_matrix_alloc(columns,columns);
    gsl_matrix *V = gsl_matrix_alloc(columns,columns);
    gsl_vector *S =gsl_vector_alloc(columns);
    gsl_vector *xx =gsl_vector_alloc(columns);  // x is in the row space, therefore it is a vector in Rcolumns, with at most row dimensions.
    gsl_vector *b =gsl_vector_alloc(rows);     // b is in the column space of A.
    gsl_vector_set_all( b,0 );
    gsl_vector *work=gsl_vector_alloc(columns);
    int rsvd;		// A must have more rows than columns or the same num
    int rsvds;               //(A,V,S,work), on output A is replaced by U  (A=USVt)
    rsvd=gsl_linalg_SV_decomp(Occupancy,V,S,work);
    if(gsl_vector_get(S,2) < .1 ) {
        cout << "gsl_matrix_get(V,2,j)   = " << gsl_matrix_get(V,2,j)  << endl;
    }
    for(int i =0; i < columns;i++){  // the number of fix_pars must be aligned with the seq2e main file..
    cout << "gsl_vector_get(S,i)  , singular values  = " << gsl_vector_get(S,i)  << endl;
}
gsl_matrix_free( Occupancy );
gsl_matrix_free( X );
gsl_matrix_free( V );
gsl_vector_free( xx );
gsl_vector_free( b );
gsl_vector_free( S );
}
int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result )   // const
{
    cout << "Start minimization simplex" << endl;
    vector< double > pars;
    printPar(par_model );
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators );
    cout << pars.size() <<  "if no num to left pars is empty:" <<endl;
    int pars_size = pars.size();
    fix_pars.clear();
    free_pars.clear();
    printPar( par_model );
    return 1;
    for( int index = 0; index < pars_size; index++ ){
        if( indicator_bool[ index ] ){
            free_pars.push_back( pars[ index ]);
        }
        else{
            fix_pars.push_back( pars[ index ] );
        }
    }
    vector < double >free_parste(0) ;
    free_parste= free_pars;
    vector < double >fix_parste(0);
    fix_parste=fix_pars;
    pars.clear();
    int free_par_counterte = 0;
    int fix_par_counterte = 0;
    for( int index = 0; index < pars_size; index ++ ){
        if( indicator_bool[ index ] ){
            pars.push_back( free_parste[ free_par_counterte ++ ]);
        }
        else{
            pars.push_back( fix_parste[ fix_par_counterte ++ ]);
        }
    }
    par_result = ExprPar( pars, coopMat, actIndicators, repIndicators );
    return 1;
    pars.clear();
    pars = free_pars;
    for (int i = 0; i < pars.size() ; i++ ) {
    }
    gsl_multimin_function my_func;
    my_func.f = &gsl_obj_f3;
    my_func.n = pars.size();
    my_func.params = (void*)this;    // params holds par_model and so does pars.., but these are essentially distinct objects, so one could change but the other not
    gsl_vector* x = vector2gsl( pars );
    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
    gsl_vector* ss = gsl_vector_alloc( my_func.n );
    gsl_vector_set_all( ss, 1.0 );
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc( T, my_func.n );
    gsl_multimin_fminimizer_set( s, &my_func, x, ss ); // here is where the problem is  101
    cout << "ssssssssssssssssssssssssssssssssssssssssssssssss" << endl;
    size_t iter = 0;
    int status;
    double size;
    printPar( par_model );
    do {
        double f_prev = iter ? s->fval : 1.0E6;     // the function starts with some very large number
        iter++;
        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_fminimizer_iterate( s );
        if ( status ) { cout<< " gsl_multimin_iterate positve: " << endl; break;}
        free_pars = gsl2vector( s-> x);
        pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < pars_size; index ++ ){
            if( indicator_bool[ index ] ){
                pars.push_back( free_pars[ free_par_counter ++ ]);
            }
            else{
                pars.push_back( fix_pars[ fix_par_counter ++ ]);
            }
        }
        ExprPar par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators );
        if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) {cout << "testpar failed: " <<endl;  break;}
        double f_curr = s->fval;            // 525 this causes the parameters to be stable or unstable if uncommented
        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_test_size( size,1e-1 ); //1e-6
        if ( status == GSL_SUCCESS ) { cout << "converged to minimum at " << iter << endl; }
        cout <<"iter" << "\t" << iter << "\t" << " f_curr " << f_curr<<endl;
    } while (  status == GSL_CONTINUE && iter < nSimplexIters );
    free_pars = gsl2vector( s-> x);
    cout << pars.size() <<  "if n8 num to left pars is empty:" <<endl;
    printPar(par_model );
    pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    for( int index = 0; index < pars_size; index ++ ){
        if( indicator_bool[ index ] ){
            pars.push_back( free_pars[ free_par_counter ++ ]);
        }
        else{
            pars.push_back( fix_pars[ fix_par_counter ++ ]);
        }
    }
    par_result = ExprPar( pars, coopMat, actIndicators, repIndicators );
    obj_result = s->fval;
    cout << "End minimization simplex" << endl;
    return 0;
}
int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result )   // const
{
    cout << "Start gradient minimiz minimization" << endl;
    printPar( par_model );
    vector< double > pars;
    printPar( par_model );
    par_model.getFreePars( pars, coopMat, actIndicators, repIndicators );
    printPar( par_model );
    cout<< "nopeeeeeeeee"<<endl;
    int pars_size = pars.size();
    fix_pars.clear();
    free_pars.clear();
    for( int index = 0; index < pars_size; index++ ){
        if( indicator_bool[ index ] ){
            free_pars.push_back( pars[ index ]);
        }
        else{
            fix_pars.push_back( pars[ index ] );
        }
    }
    pars.clear();
    pars = free_pars;
    gsl_multimin_function_fdf my_func;
    my_func.f = &gsl_obj_f3;
    my_func.df = &gsl_obj_df;
    my_func.fdf = &gsl_obj_fdf;
    my_func.n = pars.size();
    my_func.params = (void*)this;
    gsl_vector* x = vector2gsl( pars );
    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc( T, my_func.n );
    double init_step = .01, tol = .0001;
    gsl_multimin_fdfminimizer_set( s, &my_func, x, init_step, tol );
    size_t iter = 0;
    int status;
    do {
        double f_prev = iter ? s->f : 1.0E6;     // the function starts with some very large number
        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s );
        if (status ) {
            break; }
            cout << " about to printPar( par_model ) " << endl;
            printPar( par_model );
            free_pars = gsl2vector( s-> x);
            pars.clear();
            int free_par_counter = 0;
            int fix_par_counter = 0;
            for( int index = 0; index < pars_size; index ++ ){
                if( indicator_bool[ index ] ){
                    pars.push_back( free_pars[ free_par_counter ++ ]);
                }
                else{
                    pars.push_back( fix_pars[ fix_par_counter ++ ]);
                }
            }
            ExprPar par_curr = ExprPar ( pars, coopMat, actIndicators, repIndicators );
            cout << " about to printPar( par_curr ) " << endl;
            printPar( par_curr );
            if ( ExprPar::searchOption == CONSTRAINED && !testPar( par_curr ) ) { cout << "testpar failed: " <<endl;  printPar( par_curr );break;}
            double f_curr = s->f;
            double delta_f = abs( f_curr - f_prev );
            cout << " deltaf " << delta_f << " f_curr " << f_curr << " iter " << iter << endl;
            if ( objOption == SSE && delta_f < .0000001 ) {cout << "min f_sse " << endl ; break;}
            if ( objOption == CORR && delta_f < min_delta_f_Corr ) { cout << "min_delta_f_corr" << endl; break;}
            if ( objOption == CROSS_CORR && delta_f < min_delta_f_CrossCorr ) break;
            status = gsl_multimin_test_gradient( s->gradient, 5e-4 ); // was 5e-4, 81611
        } while ( status == GSL_CONTINUE && iter < nGradientIters );
        free_pars = gsl2vector( s-> x);
        pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < pars_size; index ++ ){
            if( indicator_bool[ index ] ){
                pars.push_back( free_pars[ free_par_counter ++ ]);
            }
            else{
                pars.push_back( fix_pars[ fix_par_counter ++ ]);
            }
        }
        par_result = ExprPar ( pars, coopMat, actIndicators, repIndicators );
        obj_result = s->f;
        cout << " par _result " << endl;
        printPar( par_result );
        gsl_vector_free( x );
        gsl_multimin_fdfminimizer_free( s );
        return 0;
    }
    void
    print_state (size_t iter, gsl_multifit_fdfsolver * s)
    {
        printf ("iter: %3zu x = % 15.8f % 15.8f "
        "|f(x)| = %g\n",
        iter,
        gsl_vector_get (s->x, 0),
        gsl_vector_get (s->x, 1),
        gsl_blas_dnrm2 (s->f));
    }
    int ExprPredictor::gradient_minimize2( ExprPar& par_result, double& obj_result )   // const
    {
        cout << "Start gradient minimiz minimization" << endl;
        vector< double > pars;
        par_model.getFreePars3( pars, coopMat, actIndicators, repIndicators );
        cout << "pars " << pars << endl;
        int pars_size = pars.size();
        fix_pars.clear();
        free_pars.clear();
        for( int index = 0; index < pars_size; index++ ){
            if( indicator_bool[ index ] ){
                free_pars.push_back( pars[ index ]);
            }
            else{
                fix_pars.push_back( pars[ index ] );
            }
        }
        pars.clear();
        pars = free_pars;
        const gsl_multifit_fdfsolver_type *T;
        gsl_multifit_fdfsolver *s;
        int status;
        unsigned int i, iter = 0;
        const size_t n = nSeqs()*nConds();
        const size_t p = pars.size();
        gsl_matrix* jacobian = NULL;
        gsl_vector* g = NULL;
        gsl_matrix *covar = NULL;
        cout << " n " << n << " p " << p << endl;
        gsl_multifit_function_fdf my_func;
        my_func.f = &gsl_obj_f3_fit;
        cout << " aoub " << endl;
        my_func.df = &gsl_obj_df_fit;
        cout << " aee " << endl;
        my_func.fdf =&gsl_obj_fdf_fit;
        cout << " noo " << endl;
        my_func.n = n ;
        my_func.p = p ;
        my_func.params = (void*)this;
        gsl_vector* xx = vector2gsl( pars );
        T = gsl_multifit_fdfsolver_lmder;
        s = gsl_multifit_fdfsolver_alloc (T, n, p);
        cout << " ok a " << endl;
        gsl_multifit_fdfsolver_set (s, &my_func, xx);
        const gsl_matrix* a;//(s->J);
        cout << " rint " << endl;
        bool aa=1;
        do
        {
            iter++;
            cout << "here5 " << endl;
            status = gsl_multifit_fdfsolver_iterate(s);
            printf ("status = %s\n", gsl_strerror (status));
            jacobian = gsl_matrix_alloc(n, p);
            gsl_multifit_fdfsolver_jac(s, jacobian);
            g = gsl_vector_alloc( jacobian->size2 );
            const gsl_matrix* a(jacobian);
            const gsl_vector* fa(s->f);
            cout << endl;
            cout << " okinside grad minimize 3 " <<endl;
            int aaa=0;
            aaa =gsl_multifit_gradient(a,fa, g) ;
            gsl_vector_fprintf (stdout, g, "%g");
            cout << " grad above " << endl;
            double epsabs=.001;
            status= gsl_multifit_test_gradient (g, epsabs);
            status = gsl_multifit_test_delta (s->dx, s->x,1e-4, 1e-4);
            gsl_vector_fprintf (stdout, s->dx, "%g");
            gsl_vector_fprintf (stdout, s->x, "%g");
            cout << "here5b " << endl;
        }
        while (status == GSL_CONTINUE && iter < 15);
        free_pars = gsl2vector( s-> x);
        pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < pars_size; index ++ ){
            if( indicator_bool[ index ] ){
                pars.push_back( free_pars[ free_par_counter ++ ]);
            }
            else{
                pars.push_back( fix_pars[ fix_par_counter ++ ]);
            }
        }
        par_result = ExprPar ( pars, coopMat, actIndicators, repIndicators,aa );
        covar = gsl_matrix_alloc (p, p);
        gsl_multifit_covar (jacobian, 0.0, covar);
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
        {
            double chi = gsl_blas_dnrm2(s->f);
            double dof = n - p;
            double c = GSL_MAX_DBL(1, chi / sqrt(dof));
            obj_result = chi;
            printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
            printf("chisq = %g\n", chi);
            printf ("par1 = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        }
        printf ("status = %s\n", gsl_strerror (status));
        cout << "covar = JtJ inverse, sqrt((JtJ)^-1) " << endl;   /// why does this say inverse?
        for( int ii =0 ; ii < covar->size1; ii++) {
            for( int j =0 ; j < covar->size2; j++) {
                cout <<gsl_matrix_get(covar,ii,j) << '\t';
            }
            cout << endl;
        }
        gsl_multifit_fdfsolver_free (s);
        gsl_vector_free( xx );
        gsl_matrix_free( jacobian );
        gsl_matrix_free( covar );
        gsl_vector_free( g );
        return 0;
    }
    double gsl_obj_f( const gsl_vector* v, void* params )//place compmatocc in here since we cant access v elsewhere
    {
        ExprPredictor* predictor = (ExprPredictor*)params;
        vector <double> temp_free_pars = gsl2vector(v);
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
            if( predictor ->  indicator_bool[ index ]  ){
                all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );// this constructor uses two different ExprPar objects to initialize itself: all_pars that contain free_pars comes from gsl x (i.e. v) while the rest comes from params field is the gsl struct which is type cast back to a ExprPredictor object.
        if( (*predictor).clasvar == 0 ) {
            double tempobj2 = predictor->objFuncborder2( par ); return tempobj2;// 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.
        }
        else {
            double obj = predictor->objFunc2( par );  return obj;
        }
        double tempobj2 = predictor->objFuncborder2( par ); // 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.  // 111611 changed objFuncborder2 to objFuncborder, and changed the negative sign in front of compavgcorrborder.
        return tempobj2;
    }
    double gsl_obj_f3( const gsl_vector* v, void* params )//place compmatocc in here since we cant access v elsewhere
    {
        ExprPredictor* predictor = (ExprPredictor*)params;
        vector <double> temp_free_pars = gsl2vector(v);
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
            if( predictor ->  indicator_bool[ index ]  ){
                all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );// this constructor uses two different ExprPar objects to initialize itself: all_pars that contain free_pars comes from gsl x (i.e. v) while the rest comes from params field is the gsl struct which is type cast back to a ExprPredictor object.
        if( (*predictor).clasvar == 1 ) {
            double tempobj3 = predictor->objFuncborder( par ); return tempobj3;// 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.
        }
        else {
            double tempobj3 = predictor->compRMSE2( par );   //objFuncborder2( par ); //
            return tempobj3;
        }
    }
    double gsl_obj_f2( const gsl_vector* v, void* params )//place compmatocc in here since we cant access v elsewhere
    {
        ExprPredictor* predictor = (ExprPredictor*)params;
        vector <double> temp_free_pars = gsl2vector(v);
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
            if( predictor ->  indicator_bool[ index ]  ){
                all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );// this constructor uses two different ExprPar objects to initialize itself: all_pars that contain free_pars comes from gsl x (i.e. v) while the rest comes from params field is the gsl struct which is type cast back to a ExprPredictor object.
        if( (*predictor).clasvar == 1 ) {
            double tempobj2 = predictor->objFuncborder2( par ); return tempobj2;// 73011 changed objfunc2 by deleting const method.  therefore placed objfunc2 before objfunc, hence objfunc will not have to changge sequences site representation also.
        }
        else {
            double obj = predictor->objFunc( par );  return obj;
        }
    }
    void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
    {
        double step = 1.0E-3;
        numeric_deriv( grad, gsl_obj_f3, v, params, step );
    }
    void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
    {
        *result = gsl_obj_f3( v, params );
        gsl_obj_df( v, params, grad );
    }
    int gsl_obj_f3_fit( const gsl_vector* v, void* params, gsl_vector * f) //place compmatocc in here since we cant access v elsewhere
    {
        ExprPredictor* predictor = (ExprPredictor*)params;
        vector <double> temp_free_pars = gsl2vector(v);
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        vector<double> fix_pars(predictor->getfixpars());
        for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
            if( predictor ->  indicator_bool[ index ]  ){
                all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
            }
        }
        ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() );
        vector<double>ff(0);
        predictor -> compf( par, ff);
        for ( int i = 0; i < ff.size(); i++ ) gsl_vector_set( f, i, ff[ i ] );
        return GSL_SUCCESS;
    }
    int gsl_obj_df_fit( const gsl_vector* v, void* params, gsl_matrix * J )
    {
        ExprPredictor* predictor = (ExprPredictor*)params;
        vector <double> temp_free_pars = gsl2vector(v);
        vector < double > all_pars;
        all_pars.clear();
        int free_par_counter = 0;
        int fix_par_counter = 0;
        vector< double > fix_pars(predictor->getfixpars());
        for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
            if( predictor ->  indicator_bool[ index ]  ){
                all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
            }
            else{
                all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
            }
        }
        bool a = 1;
        ExprPar par22( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(),a );
        vector< double > pars;
        pars=all_pars;
        pars.clear();
        pars =  temp_free_pars;
        double step = .2;
        vector< double > initilizer( pars.size() );
        vector< vector< double > > init2( pars.size(), initilizer );
        int count = 0; // used to indicate data point within bords*seq
        for(int i=0; i < predictor->nSeqs(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
        for(int j=0; j < predictor->nConds(); j++ ) {
            for (int jj = 0; jj < pars.size() ; jj++ ) {
                vector< double > parsjuh =pars;
                vector< double > parsjdh= pars;
                parsjuh[jj] = pars[jj] + step;
                parsjdh[jj] = pars[jj] - step;
                vector < double > all_pars;
                all_pars.clear();
                int free_par_counter = 0;
                int fix_par_counter = 0;
                vector< bool > indicator_bool(predictor->getIndicatorbool(  ));
                for( int index = 0; index < indicator_bool.size(); index ++ ){
                    if(  indicator_bool[ index ]  ){
                        all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
                    }
                    else{
                        all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                    }
                }
                const	IntMatrix coopMat( predictor->getcoopmat());
                const vector< bool > actIndicators(predictor->getactIndicators());
                const vector< bool > repIndicators(predictor->getrepIndicators());
                ExprPar parjuh = ExprPar( all_pars, coopMat, actIndicators, repIndicators , a);
                all_pars.clear();
                free_par_counter = 0;
                fix_par_counter = 0;
                for( int index = 0; index < indicator_bool.size(); index ++ ){
                    if(  indicator_bool[ index ]  ){
                        all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
                    }
                    else{
                        all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                    }
                }
                ExprPar parjdh = ExprPar( all_pars, coopMat, actIndicators, repIndicators, a );
                vector< double > fOcc(0);
                count = predictor->nSeqs();
                double pjuh = predictor->compRMSE4(parjuh, i,j); // func->predictExpr3(seqSitesm1[i], 5, concs );
                double pjdh =  predictor->compRMSE4(parjdh, i,j); //func->predictExpr3(seqSitesm1[i], 5, concs );
                gsl_matrix_set( J,predictor->nConds()*i+j, jj, ( pjuh- pjdh )/(2.0*step) );  // i is seqs, jj is parameters
                parsjuh.clear();
                parsjdh.clear();
            } //for j nconds
        }  // for jj npars
    }// for i nseqs
    return GSL_SUCCESS;
}
int gsl_obj_fdf_fit(  const gsl_vector* v, void* params, gsl_vector* f, gsl_matrix * J)
{
    ExprPredictor* predictor = (ExprPredictor*)params;
    vector <double> temp_free_pars = gsl2vector(v);
    vector < double > all_pars;
    all_pars.clear();
    int free_par_counter = 0;
    int fix_par_counter = 0;
    vector<double> fix_pars(predictor->getfixpars());
    for( int index = 0; index < predictor -> indicator_bool.size(); index ++ ){
        if( predictor ->  indicator_bool[ index ]  ){
            all_pars.push_back( temp_free_pars[ free_par_counter ++ ]  );
        }
        else{
            all_pars.push_back( predictor ->  fix_pars[ fix_par_counter ++ ]  );
        }
    }
    bool a=1;
    ExprPar par( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators() ,a);
    cout << " inside here " << endl;
    vector<double >ff(0);
    predictor -> compf( par, ff);
    for ( int i = 0; i < ff.size(); i++ ) gsl_vector_set( f, i, ff[ i ] );
    vector< double > pars;
    pars=all_pars;
    pars.clear();
    pars =  temp_free_pars;
    double step = .2;
    vector< double > initilizer( pars.size() );
    vector< vector< double > > init2( pars.size(), initilizer );
    int count = 0; // used to indicate data point within bords*seq
    for(int i=0; i < predictor->nSeqs(); i++ ) {   // this loops over sequences, i.e. , i represents a sequence
    for(int j=0; j < predictor->nConds(); j++ ) {
        for (int jj = 0; jj < pars.size() ; jj++ ) {
            vector< double > parsjuh =pars;
            vector< double > parsjdh= pars;
            parsjuh[jj] = pars[jj] + step;
            parsjdh[jj] = pars[jj] - step;
            vector < double > all_pars;
            all_pars.clear();
            int free_par_counter = 0;
            int fix_par_counter = 0;
            vector< bool > indicator_bool(predictor->getIndicatorbool(  ));
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjuh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            const	IntMatrix coopMat( predictor->getcoopmat());
            const vector< bool > actIndicators(predictor->getactIndicators());
            const vector< bool > repIndicators(predictor->getrepIndicators());
            ExprPar parjuh =ExprPar( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(),a );
            all_pars.clear();
            free_par_counter = 0;
            fix_par_counter = 0;
            for( int index = 0; index < indicator_bool.size(); index ++ ){
                if(  indicator_bool[ index ]  ){
                    all_pars.push_back( parsjdh[ free_par_counter ++ ]  );
                }
                else{
                    all_pars.push_back( fix_pars[ fix_par_counter ++ ]  );
                }
            }
            ExprPar parjdh =ExprPar( all_pars, predictor->getCoopMat(), predictor->getActIndicators(), predictor->getRepIndicators(),a );
            vector< double > fOcc(0);
            count = predictor->nSeqs();
            double pjuh = predictor->compRMSE4(parjuh, i,j); // func->predictExpr3(seqSitesm1[i], 5, concs );
            double pjdh =  predictor->compRMSE4(parjdh, i,j); //func->predictExpr3(seqSitesm1[i], 5, concs );
            gsl_matrix_set( J,predictor->nConds()*i+j, jj, ( pjuh- pjdh )/(2.0*step) );  // i is seqs, jj is parameters
            parsjuh.clear();
            parsjdh.clear();
        } //for j nconds
    }  // for jj npars
}// for i nseqs
for ( int i = 0; i < predictor->nSeqs(); i++ ) {
    for ( int j = 0; j < predictor->nConds(); j++ ) {
    }
}
return GSL_SUCCESS;
}
