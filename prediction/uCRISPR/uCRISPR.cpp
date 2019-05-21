#include <iostream>
using namespace std;
#include <fstream> //printf()
#include <cstdlib> //atoi(), atof(), exit()
#include <iomanip> //setprecision()
#include <cmath>   //exp()
#include <ctime>
#include <sstream> //string
#include <vector>
#include <map>

string DATA="/home/yjf/off-target/off-target/tool/uCRISPR/data/";   //absolute path
string RNAstructure="/home/yjf/off-target/off-target/tool/uCRISPR/RNAstructure/exe/EnsembleEnergy"; //absolute path

vector<string> string_to_vector_string(string sline, int outlet) {
   istringstream ss (sline);
   string buf;
   vector<string> token;
   while(ss >> buf) token.push_back(buf);
   
   if(outlet == 1) {
      for(vector<string>::iterator lt = token.begin(); lt != token.end(); lt ++)
         cout << *lt << " ";
      cout << endl;
   }
   
   return token;
}

map<string,double> parameters_on;
map<string,double> parameters_off;
int ReadParameters() {
   string sline;
   vector<string> stoken;

   string fstr=DATA+"parameters-ontarget.dat";
   ifstream inP (fstr.c_str());
   if(!inP.is_open()) {
     cout<<"Can't open file "<<fstr<<endl;
     return 1;
   }
   while(getline(inP,sline)) {
     stoken=string_to_vector_string(sline,0);
     if(stoken.size()!=2) continue;
     string name=stoken[1];
     double value=atof(stoken[0].c_str());
     parameters_on.insert(pair<string,double>(name,value));
   }
   inP.close(); inP.clear();

   fstr=DATA+"parameters-offtarget.dat";
   inP.open(fstr.c_str());
   if(!inP.is_open()) {
     cout<<"Can't open file "<<fstr<<endl;
     return 1;
   }
   while(getline(inP,sline)) {
     stoken=string_to_vector_string(sline,0);
     if(stoken.size()!=2) continue;
     string name=stoken[1];
     double value=atof(stoken[0].c_str());
     parameters_off.insert(pair<string,double>(name,value));
   }
   inP.close(); inP.clear();

   return 0;
}

string CheckSequence(string sequence) {
   string seq="";
   for(int ii=0;ii<sequence.size();ii++) {
     char nt=sequence[ii];
     if(nt=='A'||nt=='a') seq+="A";
     else if(nt=='G'||nt=='g') seq+="G";
     else if(nt=='C'||nt=='c') seq+="C";
     else if(nt=='U'||nt=='u') seq+="T";
     else if(nt=='T'||nt=='t') seq+="T";
     else 
       return "-1";
   }

   return seq;
}

map<string,double> EsgRNA_reserved;
double GetSgRNAEnergy(string sequence) {
   string sline;
   string fstr=DATA+"sgRNA.fasta";
   ifstream inS (fstr.c_str());
   getline(inS,sline);
   getline(inS,sline);
   string sgRNA=sline.substr(20);
   string gRNA=sequence;
   while(gRNA.find("T")!=string::npos) {
     int found=gRNA.find("T");
     gRNA.replace(found,1,"U");
   }
   sgRNA=gRNA+sgRNA;

   if(EsgRNA_reserved.count(sgRNA))
     return EsgRNA_reserved[sgRNA];

   ofstream outS ("gRNA.seq");
   outS<<sgRNA<<endl;
   outS.close(); outS.clear();
   string cmd=RNAstructure+" gRNA.seq -s --sequence > gRNA.eng";
   system(cmd.c_str());
   ifstream inE ("gRNA.eng");
   getline(inE,sline);
   inE.close(); inE.clear();
   int found=sline.find(":");
   string str=sline.substr(found+1);
   vector<string> stoken=string_to_vector_string(str,0);
   double energy=atof(stoken[0].c_str());

   cmd="rm gRNA.seq gRNA.eng";
   system(cmd.c_str());  

   EsgRNA_reserved.insert(pair<string,double>(sgRNA,energy));

   return energy;   
}

double CalculateOntargetScore(string sequence) {
   string dna="";
   string rna="";
   string pam="";
   if(sequence.size()==23) {
     dna=sequence;
     rna=sequence.substr(0,20);
     pam=sequence.substr(21);
   }  
   else if(sequence.size()==30) {
     dna=sequence.substr(0,25);
     dna+=sequence.substr(27);
     rna=sequence.substr(4,20);
     pam=sequence.substr(25,2);
   }
   
   double E_C=0.0;
   if(sequence.size()==23) {
     for(int ii=0;ii<20;ii++) {
       string pname=dna.substr(ii,2)+"_"+to_string(ii+1);
       E_C+=parameters_on[pname];
     }
   }
   if(sequence.size()==30) {
     for(int ii=-3;ii<24;ii++) {
       string pname=dna.substr(ii+3,2)+"_"+to_string(ii);
       E_C+=parameters_on[pname];
     }
   }

   double EsgRNA=GetSgRNAEnergy(rna);
   EsgRNA*=parameters_on["wdG"];

   double Epam=1.0;
   if(pam=="GG") Epam=1.0;
   else if(pam=="AG") Epam=0.259259259;
   else if(pam=="CG") Epam=0.107142857;
   else if(pam=="GA") Epam=0.069444444;
   else if(pam=="GC") Epam=0.022222222;
   else if(pam=="GT") Epam=0.016129032;
   else if(pam=="TG") Epam=0.038961039;
   else Epam=0.0;

   double score=exp(E_C+EsgRNA)*Epam;

   return score;
}

double CalculateOfftargetScore(string seq_wt, string seq_ot) {
   string rna=seq_wt;
   string dna=seq_ot; 
   double Emm=0.0;
   string smatch="";
   for(int kk=0;kk<20;kk++) {
     if(rna[kk]!=dna[kk]) {
       string pname=rna.substr(kk,1)+dna.substr(kk,1)+"_"+to_string(kk+1);
       Emm+=parameters_off[pname];
       smatch+="0";
     }
     else
       smatch+="1";
   }
   vector<string> cmm;
   smatch+="1";
   while(1) {
     if(smatch.find("0")==string::npos) break;
     int found0=smatch.find_first_of("0");
     if(found0!=0)
       smatch=smatch.substr(found0);
     int found1=smatch.find_first_of("1");
     string str=smatch.substr(0,found1);
     if(str.size()>1) 
       cmm.push_back(str);
     smatch=smatch.substr(found1);
   }
   for(int kk=0;kk<cmm.size();kk++) {
     int len=cmm[kk].size();
     if(len>5) len=5;
     string pname="CMM_"+to_string(len);
     Emm+=parameters_off[pname];
   }
      
   double E_C=0.0;
   for(int kk=0;kk<19;kk++) {
     string dirna=rna.substr(kk,2)+"_"+to_string(kk+1);
     string didna=dna.substr(kk,2)+"_"+to_string(kk+1);
     if(dirna==didna)
       E_C+=parameters_on[didna];
     else if(dirna[0]==didna[0]||dirna[1]==didna[1]) {
       double wg1=parameters_off["wg1"];
       double wt1=parameters_off["wt1"];
       Emm+=wg1*parameters_on[dirna]+wt1*parameters_on[didna];
     }
     else {
       double wg2=parameters_off["wg2"];
       double wt2=parameters_off["wt2"];
       Emm+=wg2*parameters_on[dirna]+wt2*parameters_on[didna];
     }
   }
   string pname=dna.substr(19,2)+"_20";
   E_C+=parameters_on[pname];

   double EsgRNA=GetSgRNAEnergy(rna);
   EsgRNA*=parameters_on["wdG"];

   string pam=dna.substr(21);
   double Epam=1.0;
   if(pam=="GG") Epam=1.0;
   else if(pam=="AG") Epam=0.259259259;
   else if(pam=="CG") Epam=0.107142857;
   else if(pam=="GA") Epam=0.069444444;
   else if(pam=="GC") Epam=0.022222222;
   else if(pam=="GT") Epam=0.016129032;
   else if(pam=="TG") Epam=0.038961039;
   else Epam=0.0;

   double score=exp(E_C+Emm+EsgRNA)*Epam;

   return score;
}

int main (int argc,char **argv) {
   int hcheck=0;
   if(argc!=3) hcheck=1;
   else {
     string cline=argv[1];
     if(cline!="-on"&&cline!="-off") hcheck=1;
   }
   if(hcheck) {
     cout<<"Usage: uCRISPR [options]"<<endl;
     cout<<"Options:"<<endl;
     cout<<"   -h        #Display help information"<<endl;
     cout<<"   -on file  #Evaluate on-target avtivities for on-target sites in file,"<<endl;
     cout<<"              one site (23-mer or 30-mer sequence) per line."<<endl;
     cout<<"              See \"./example/OntargetSite.dat\" for an example."<<endl;
     cout<<"   -off file #Evaluate off-target efficiencies for off-target sites in file,"<<endl;
     cout<<"              one site (20-mer/23-mer wild type sequence and 23-mer off target"<<endl; 
     cout<<"              sequence) per line in file. See \"./example/OfftargetSite.dat\" for an example."<<endl;
     cout<<"=============================================================================================="<<endl;
     cout<<"Examples:"<<endl;
     cout<<"uCRISPR -on  ./example/OntargetSite.dat   #Evaluate on-target sites given in file \"OntargetSite.dat\","<<endl;
     cout<<"                                           predicted scores are shown on the screen."<<endl;
     cout<<"uCRISPR -off ./example/OfftargetSite.dat  #Evaluate off-target sites given in file \"OfftargetSite.dat\","<<endl;
     cout<<"                                           predicted scores are shown on the screen."<<endl;
     
     return 1;
   }
   
   string flag=argv[1];
   string filename=argv[2];

   string sline;
   vector<string> Vseq;
   ifstream inD (filename.c_str());
   if(!inD.is_open()) {
     cout<<"Can't open file "<<filename<<endl;
     return 1;
   }
   while(getline(inD,sline)) {
     Vseq.push_back(sline);
   }
   inD.close(); inD.clear();

   if(ReadParameters()) return 1;

   vector<string> stoken;
   vector<string> Vout;
   if(flag=="-on") {
     Vout.push_back("#OntargetSite uCRISPR_score");
     for(int ii=0;ii<Vseq.size();ii++) {
       if(Vseq[ii].find("#")!=string::npos) {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Found symbol '#'."<<endl;
         continue;
       }
       stoken=string_to_vector_string(Vseq[ii],0);
       string sequence=stoken[0];
       string seq_on=CheckSequence(sequence);
       if(seq_on=="-1") {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Unknown nucleotide identity";
         cout<<"(Only A/C/G/T/U/a/c/g/t/u accepted)"<<endl;
         continue;
       }
       if(seq_on.size()!=23&&seq_on.size()!=30) {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Ontarget site should be 23-mer or 30-mer long.";
         continue;
       }
       double score=CalculateOntargetScore(seq_on);
       string str=sequence+" "+to_string(score);
       Vout.push_back(str);
     }
   }
   else if(flag=="-off") {
     Vout.push_back("#Wildtype_sequence Offtarget_sequence uCRISPR_score");
     for(int ii=0;ii<Vseq.size();ii++) {
       if(Vseq[ii].find("#")!=string::npos) {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Found symbol '#'."<<endl;
         continue;
       }
       stoken=string_to_vector_string(Vseq[ii],0);
       string sequence_0=stoken[0];
       if(sequence_0.size()!=20&&sequence_0.size()!=23) {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Wildtype sequence should be 20-mer or 23-mer long.";
         continue;
       }
       sequence_0=sequence_0.substr(0,20);
       string seq_wt=CheckSequence(sequence_0);
       if(seq_wt=="-1") {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Unknown nucleotide identity";
         cout<<"(Only A/C/G/T/U/a/c/g/t/u accepted)"<<endl;
         continue;
       }
       string sequence_1=stoken[1];
       if(sequence_1.size()!=23) {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Offtarget sequence should be 23-mer long.";
         continue;
       }
       string seq_ot=CheckSequence(sequence_1);
       if(seq_ot=="-1") {
         cout<<"Line "<<ii+1<<" \""<<Vseq[ii]<<"\" is skipped: Unknown nucleotide identity";
         cout<<"(Only A/C/G/T/U/a/c/g/t/u accepted)"<<endl;
         continue;
       }
       double score=CalculateOfftargetScore(seq_wt,seq_ot);
       string str=sequence_0+" "+sequence_1+" "+to_string(score);
       Vout.push_back(str);
     }
   }

   for(int ii=0;ii<Vout.size();ii++)
     cout<<Vout[ii]<<endl;

   return 0;
}