
//#include <time.h>
# include <stdio.h>
#include <iostream>
//#include <string>
# include <math.h>
#include <vector>
using namespace std;

const int NumSensor=10;	// totale sensori

///////////////////////////////////////////////////////////////////
//Utilita'
double Segno(double numero);
double TrovaMinimo(double *Array, int MaxDim);
int TrovaIndiceDelMinimo(double *Array, int MaxDim);
int TrovaIndiceDelMassimo(double *Array, int MaxDim);
double TrovaMassimo(double *Array, int MaxDim);



double NumeroRand(double max);





int ClassificatoreLineare(double x, double y);
int ClassificatoreLineare3D(double x, double y, double z);
int ClassificatoreLineare6D(double x, double y, double z, double x2, double y2, double z2);
int ClassificatoreLineare12D(double x, double y, double z, double x2, double y2, double z2, 
												double x3, double y3, double z3, double x4, double y4, double z4, int classificatore=0);
int ClassificatoreQ(double m, double v);														
int ClassificatoreRettangolo(double Dt_mmin, double A_mmin, double Dt_mMax, double A_mMax, double Dt_mSample, double A_mSample);
int ClassificatoreRettangolo3D(double Dt_mmin, double A_mmin, double Q_mmin, double Dt_mMax, double A_mMax, double Q_mMax, 
							   double Dt_mSample,  double A_mSample,  double Q_mSample);
int ClassificatoreRettangolo6D(double Dt_mmin, double A_mmin, double Q_mmin, double Dt_mMax, double A_mMax, double Q_mMax,
							   double Dt_vmin, double A_vmin, double Q_vmin, double Dt_vMax, double A_vMax, double Q_vMax,
							   double Dt_mSample,  double A_mSample,  double Q_mSample,
							   double Dt_vSample,  double A_vSample,  double Q_vSample);
int ClassificatoreRettangolo12D(double Dt_mmin, double A_mmin, double Q_mmin, double Dt_mMax, double A_mMax, double Q_mMax,
							   double Dt_vmin, double A_vmin, double Q_vmin, double Dt_vMax, double A_vMax, double Q_vMax,
							   double Dt_smin, double A_smin, double Q_smin, double Dt_sMax, double A_sMax, double Q_sMax,
							   double Dt_kmin, double A_kmin, double Q_kmin, double Dt_kMax, double A_kMax, double Q_kMax,
							   double Dt_mSample,  double A_mSample,  double Q_mSample,
							   double Dt_vSample,  double A_vSample,  double Q_vSample,
							   double Dt_sSample,  double A_sSample,  double Q_sSample,
							   double Dt_kSample,  double A_kSample,  double Q_kSample);			
void Permuta(int volte, int Dim, double *x, double *y, double *z);
double ErroreLineare12D(int Size, int classificatore=0, double *perf=NULL);
double ErroreRules(int volte, int div_anomaly, double *perf=NULL);
double ErroreLineare12D_ripetuto_amaggioranza(int Size, int classificatore, int num_campioni);

double ErroreRules_ripetuto_amaggioranza(int Size, int num_campioni);
int VotingLineare12D(double x, double y, double z, double x2, double y2, double z2, 
												double x3, double y3, double z3, double x4, double y4, double z4, double *pesi);
double ErroreVoting12DD(int Size, double *pesi);
void Voting12DCFSQP(double *pesi);
void obj32_Voting12D(int nparam, int j, double *a_csfqp,double *fj, void *cd);
void cntr32_Voting12D(int nparam,int j,double *a_csfqp,double *gj,void *cd);

// copia interna del db DtAQSK_stat.txt 
double *Dt_m; double *A_m; double *Q_m; 
double *Dt_v; double *A_v; double *Q_v;
double *Dt_s; double *A_s; double *Q_s;
double *Dt_k; double *A_k; double *Q_k; 
int *DtAQSK_out; 
int indiceDtAQSK;

double mapminmax(double x, double xmin, double xMax);
double mapminmax_reverse(double y, double xmin, double xMax);
double ymin=-1, ymax=1.0;
double Dt_mMax2, Dt_mmin2, A_mMax2, A_mmin2, Q_mMax2, Q_mmin2,
           Dt_vMax2, Dt_vmin2, A_vMax2, A_vmin2, Q_vMax2, Q_vmin2,
		   Dt_sMax2, Dt_smin2, A_sMax2, A_smin2, Q_sMax2, Q_smin2,
		   Dt_kMax2, Dt_kmin2, A_kMax2, A_kmin2, Q_kMax2, Q_kmin2;

int K=1e4*3; //1e5*3; //K is the size of the training set (divided 50% between db tunnel and db legitimate (i.e., without tunnel))
int n_s=1e3*1;//1e3*1.2; //n_s is the # of rows to build the statistics of the feature vector
bool attivaMix=1;
int contatoreMix=1; //1 50% - 9 10% - 99 1% - 999 0.1% - 5 20%

FILE *Out2; int PPP=0, NNN=0, div_anomaly; 

////////////////////////////////////////////////////////////////////
void main(){	
		int i, j, p, cerca; srand(90); 
		
	const int Dimdns_area2=267433;
	const int Dimdns_ieiit=158725;
	const int Dimdns_ieiit2=141338;
	const int Dimarea2ieiit=426157;	
	
	const int Dimdns2tcpwget=181822;
	const int Dimdns2tcpp2p=184977;
	const int Dimdns2tcpssh=174991;
	const int Dimdns2tcpall=541788;

	const int Dimiodine_wget=102538;
	const int Dimiodine_p2p=105981;
	const int Dimiodine_ssh=99970;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	const int dimlog= Dimdns_ieiit-1 ;	// #righe file-1 
	//const int dimlog= Dimdns_area2-1 ;	
	const int dimlog2= Dimdns2tcpp2p-1 ;	// #righe file-1 
	//const int dimlog2= Dimiodine_p2p-1 ;	// #righe file-1 
	
	FILE *Out;
	FILE *dns=fopen("dns_ieiit.txt","r"); // dns_ieiit2			
	//FILE *dns=fopen("dns_area2.txt","r");	
	FILE *dns2=fopen("dns2tcpp2p.txt","r");	// in case of ssh, replace dns2tcpp2p with dns2tcpssh  	
	//iodine_p2p iodine_ssh iodine_wget
	//FILE *dns2=fopen("iodine_ssh.txt","r");		
		
	div_anomaly=1; 
	
	double t1,t2;

	//const int elementiperriga=6;
	double *sizeQ;
	double *sizeA;
	double *tQ;
	double *tA;
	int *N_filelog;
	int *IdDns;			
	char **src;
	char **dst;
	double *Dt;

	double *sizeQ2;
	double *sizeA2;
	double *tQ2;
	double *tA2;
	int *N_filelog2;
	int *IdDns2;			
	char **src2;
	char **dst2;
	double *Dt2;

	double *DtMix;
	double *AMix;
	double *QMix;

	char local[20];

	sizeQ=new double [dimlog];	
	sizeA=new double [dimlog];	
	tQ=new double [dimlog];	
	tA=new double [dimlog];	
	N_filelog=new int [dimlog];	
	IdDns=new int [dimlog];	
	Dt=new double [dimlog];	

	sizeQ2=new double [dimlog2];	
	sizeA2=new double [dimlog2];	
	tQ2=new double [dimlog2];	
	tA2=new double [dimlog2];	
	N_filelog2=new int [dimlog2];	
	IdDns2=new int [dimlog2];	
	Dt2=new double [dimlog2];	

	DtMix=new double [dimlog+dimlog2];	
	AMix=new double [dimlog+dimlog2];
	QMix=new double [dimlog+dimlog2];

	src=new char *[dimlog];
	for(p=0;p<dimlog;p++)
		src[p]=new char[20];
	dst=new char *[dimlog];
	for(p=0;p<dimlog;p++)
		dst[p]=new char[20];	

	src2=new char *[dimlog2];
	for(p=0;p<dimlog2;p++)
		src2[p]=new char[20];
	dst2=new char *[dimlog2];
	for(p=0;p<dimlog2;p++)
		dst2[p]=new char[20];	

	////////////////////////////////////////// varie sperimentazioni di lettura stringhe //////////////////////////////////////////
	/*char   string[100];
	int p;
	fscanf(dns,"%d\t", &p);
	fscanf(dns,"%s", &string);
	*/
	/*
	string s2 ( "ma vaffambagno" );   
	printf("\n\n%s",s2.data());
	printf("\n\n"); return;
	string src[dimlog]; //serve #include <string>
	string dst[dimlog];
	*/
	////////////////////////////////////////// END varie sperimentazioni di lettura stringhe //////////////////////////////////////////
	
	//Salto la prima riga... dns
	char local2[20];
	for(i=0;i<7;i++)
		fscanf(dns,"%s\t", &local2);
	fscanf(dns,"%s\n", &local2);

	//Salto la prima riga... dns2
	for(i=0;i<7;i++)
		fscanf(dns2,"%s\t", &local2);
	fscanf(dns2,"%s\n", &local2);

	for(i=0;i<dimlog;i++){

		fscanf(dns,"%d\t", &N_filelog[i]);
		fscanf(dns,"%d\t", &IdDns[i]);
		fscanf(dns,"%lf\t", &sizeQ[i]);
		fscanf(dns,"%lf\t", &tQ[i]); if(i==0) t1=tQ[i]; if(i==(dimlog-1)) t2=tQ[i]; //piglio alla fine il tempo della Q perchè a volte la A è a tempo 0.0
		fscanf(dns,"%lf\t", &sizeA[i]);
		fscanf(dns,"%lf\t", &tA[i]); 			

		fscanf(dns,"%s\t", &local);
		//printf("\n%s",local);		
		strncpy (src[i],local,20);
		//printf("\n%s",src[i]);		

		fscanf(dns,"%s\n", &local);
		strncpy (dst[i],local,20);		
	}//END for su dimlog
	printf("\n%.2f\t%2f\t",t1,t2); 
	printf("\nTempo di monitoraggio=%.2f\tFrequenza=%2f",(t2-t1),dimlog/(t2-t1)); 
	printf("\nTempo per 1e3 campioni=%.2f [s]",((t2-t1)*1e3)/dimlog); 
	printf("\nTempo per 1e4 campioni=%.2f [m]",(((t2-t1)*1e4)/dimlog)/60); 
	printf("\nTempo per 1e5 campioni=%.2f [m]",(((t2-t1)*1e5)/dimlog)/60); 
	printf("\n\n");	//system("pause");

	//printf("\n----------------------- * -----------------------");
	for( i=dimlog-5; i<0; i++ ){//stampa della lettura
		
		printf("\n");		
		printf("%d\t",N_filelog[i]);
		printf("%d\t",IdDns[i]);
		printf("%f\t",sizeQ[i]);
		printf("%f\t",tQ[i]);
		printf("%f\t",sizeA[i]);
		printf("%f\t",tA[i]);
		printf("%s\t",src[i]);
		printf("%s\t",dst[i]);

	}//END for su dimlog di stampa della lettura 
	//system("pause");

	for(i=0;i<dimlog2;i++){

		fscanf(dns2,"%d\t", &N_filelog2[i]);
		fscanf(dns2,"%d\t", &IdDns2[i]);
		fscanf(dns2,"%lf\t", &sizeQ2[i]);
		fscanf(dns2,"%lf\t", &tQ2[i]); if(i==0) t1=tQ2[i]; if(i==(dimlog2-1)) t2=tQ2[i]; //piglio alla fine il tempo della Q perchè a volte la A è a tempo 0.0		
		fscanf(dns2,"%lf\t", &sizeA2[i]);
		fscanf(dns2,"%lf\t", &tA2[i]); 		

		fscanf(dns2,"%s\t", &local);
		//printf("\n%s",local);		
		strncpy (src2[i],local,20);
		//printf("\n%s",src[i]);		

		fscanf(dns2,"%s\n", &local);
		strncpy (dst2[i],local,20);		
	}//END for su dimlog	
	printf("\n%.2f\t%2f\t",t1,t2); 
	printf("\nTempo di monitoraggio=%.2f\tFrequenza=%2f",(t2-t1),dimlog2/(t2-t1)); 
	printf("\nTempo per 1e3 campioni=%.2f [s]",((t2-t1)*1e3)/dimlog2); 
	printf("\nTempo per 1e4 campioni=%.2f [m]",(((t2-t1)*1e4)/dimlog2)/60); 
	printf("\nTempo per 1e5 campioni=%.2f [m]",(((t2-t1)*1e5)/dimlog2)/60); 
	printf("\n\n");	//system("pause");
	
	//printf("\n----------------------- * -----------------------");
	for( i=dimlog2-5; i<0; i++ ){//stampa della lettura
		
		printf("\n");		
		printf("%d\t",N_filelog2[i]);
		printf("%d\t",IdDns2[i]);
		printf("%f\t",sizeQ2[i]);
		printf("%f\t",tQ2[i]);
		printf("%f\t",sizeA2[i]);
		printf("%f\t",tA2[i]);
		printf("%s\t",src2[i]);
		printf("%s\t",dst2[i]);

	}//END for su dimlog di stampa della lettura 
	//system("pause");

	int k=-1;
	for(i=0;i<dimlog;i++){
		if(sizeA[i]>0.0){
			k++;
			Dt[k]=tA[i]-tQ[i];
			sizeQ[k]=sizeQ[i];
			sizeA[k]=sizeA[i];			
		}
	}
	int dimlog_ripulito=k;
	printf("\n\ndimlog_ripulito=%d",dimlog_ripulito);		

	char *stringadacercare2="150.145.0.109";
	char *stringadacercare3="150.145.5.118";	

	k=-1;
	for(i=0;i<dimlog2;i++){			
		if(sizeA2[i]>0.0 
				/*&& (strncmp(src2[i],stringadacercare2,(unsigned)strlen(stringadacercare2)) != 0) && (strncmp(src2[i],stringadacercare3,(unsigned)strlen(stringadacercare3)) != 0)*/ ){
			k++;
			Dt2[k]=tA2[i]-tQ2[i];
			sizeQ2[k]=sizeQ2[i];
			sizeA2[k]=sizeA2[i];		
			//IdDns2[k]=IdDns2[i];
			/*if ( strncmp(src2[i],stringadacercare2,(unsigned)strlen(stringadacercare2)) == 0){
				printf ("\nbeccato %d",IdDns2[k]);system("pause");
			}*/	
			/*if(Dt2[k]>18)
				printf ("\n%s\t%s",src2[i],dst2[i]);*/
		}
	}
	int dimlog2_ripulito=k;
	printf("\n\ndimlog2_ripulito=%d",dimlog2_ripulito);	
	//system("pause");

	/*Out=fopen("Dt.txt","w");
	for(j=0;j<dimlog_ripulito;j++)
		fprintf(Out,"\n%lf",Dt[j]);			 
	fclose(Out);
	Out=fopen("Dt2.txt","w");
	for(j=0;j<dimlog2_ripulito;j++)
		//fprintf(Out,"\n%d\t%lf",IdDns2[j],Dt2[j]);			 
		fprintf(Out,"\n%lf",Dt2[j]);			 
	fclose(Out);*/
	//return;

	//piglio il più piccolo, devo leggere il più piccolo ed aggiungere pezzi del più grosso
	int Dim;
	if(dimlog_ripulito<=dimlog2_ripulito)
		Dim=dimlog_ripulito;
	else
		Dim=dimlog2_ripulito;

/*N.B. nel ciclo del mix devo aggiungere ai dati che prendo dal db clean i dati del tunnel quando il 'conta' raggiunge questi valori 
e generando le seguenti % di dati di tunnel dentro il db:
																				1 50% - 9 10% - 99 1%		*/
	k=-1;
	int conta=0;
	int contaMax=contatoreMix;

	int i_shift, shift=(int)NumeroRand(Dim);

	for(i=0;i<Dim;i++){		
		k++;
		if(attivaMix){ // vedi N.B. su
			DtMix[k]=Dt[i];
			AMix[k]=sizeA[i];	
			QMix[k]=sizeQ[i];	
		}else{ // senza mix devo mettere i dati del tunnel

			i_shift=i+shift-1; 
			if(i_shift>Dim-1) 
				i_shift=0;
			i_shift++;
			//i_shift=i;

			DtMix[k]=Dt2[i_shift];
			AMix[k]=sizeA2[i_shift];	
			QMix[k]=sizeQ2[i_shift];	
		}
		
		conta++;
		if(attivaMix && conta==contaMax){ 
			conta=0;
			k++;

			i_shift=i+shift-1; 
			if(i_shift>Dim-1) i_shift=0;
			i_shift++;
			//i_shift=i;

			DtMix[k]=Dt2[i_shift];
			AMix[k]=sizeA2[i_shift];	
			QMix[k]=sizeQ2[i_shift];	
		}
	}
	int dimlogMix=k;
	printf("\n\ndimlogMix=%d",dimlogMix);		

	// Dato che così ho generanto nel db mix dei dati in sequenza (e.g., 10 dell'uno 1 dell'altro, sempre in ordine), faccio un mischione randomico
	if(1) Permuta(10, dimlogMix, DtMix, AMix, QMix);

	/*Out=fopen("DtMix.txt","w");
	for(j=0;j<dimlogMix;j++)
		fprintf(Out,"\n%lf",DtMix[j]);			 
	fclose(Out);*/

	//---------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Per db e dbMix piglio a caso un punto di partenza e piglio da lì 5000 campioni e ne calcolo media e varianza, tutto ciò ripetutamente
	int NumElementi=n_s;
	int volte=K/2;	
	int posizione;
	double somma;
	double sommaquadro;
	double media;
	double var;
	int out;	
	
	/*Out=fopen("Dt_stat.txt","w");
	fprintf(Out,"m;v;g");
	fclose(Out);//per far sparire lo stesso file precedentemente creato basta anche solo aprire e chiudere senza scrivere niente

	Out=fopen("A_stat.txt","w");
	fprintf(Out,"m;v;g");
	fclose(Out);*/

	Out=fopen("DtAQSK_stat.txt","w");
	//fprintf(Out,"mDt;mA;mQ;vDt;vA;vQ;sDt;sA;sQ;kDt;kA;kQ;g");
	if(0) fprintf(Out,"mDt\tmA\tmQ\tvDt\tvA\tvQ\tsDt\tsA\tsQ\tkDt\tkA\tkQ\tg");
	//fprintf(Out,"Dt;A;Q;g");
	//fprintf(Out,"mDt;mA;mQ;vDt;vA;vQ;g");
	fclose(Out);

	if(0){
		Out=fopen("0.txt","w"); fprintf(Out,"mDt\tmA\tmQ\tvDt\tvA\tvQ\tsDt\tsA\tsQ\tkDt\tkA\tkQ"); fclose(Out); 
		Out=fopen("1.txt","w"); fprintf(Out,"mDt\tmA\tmQ\tvDt\tvA\tvQ\tsDt\tsA\tsQ\tkDt\tkA\tkQ"); fclose(Out); 
	}

	//copia interna del db Dt_stat.txt che creo ora (e passo poi al matlab per fare una prova del classificatore che piglio dal matlab)
	int dim_copiainterna=volte*2;
	double *Dt_stat_m;
	Dt_stat_m=new double [dim_copiainterna];		
	double *Dt_stat_v;
	Dt_stat_v=new double [dim_copiainterna];
	int *Dt_stat_out;
	Dt_stat_out=new int [dim_copiainterna];
	int indiceDt_stat=-1;
	//END copia interna del db Dt_stat.txt che creo ora (e passo poi al matlab per fare una prova del classificatore che piglio dal matlab)	
	
	/*Out=fopen("Dt_stat_NNInput.txt","w");
	fclose(Out);
	Out=fopen("Dt_stat_NNTarget.txt","w");
	fclose(Out);*/

	//copia interna del db DtAQSK_stat.txt 
	dim_copiainterna=volte*2;

	Dt_m=new double [dim_copiainterna];			
	A_m=new double [dim_copiainterna];	
	Q_m=new double [dim_copiainterna];
	
	Dt_v=new double [dim_copiainterna];			
	A_v=new double [dim_copiainterna];	
	Q_v=new double [dim_copiainterna];
	
	Dt_s=new double [dim_copiainterna];			
	A_s=new double [dim_copiainterna];	
	Q_s=new double [dim_copiainterna];
	
	Dt_k=new double [dim_copiainterna];			
	A_k=new double [dim_copiainterna];	
	Q_k=new double [dim_copiainterna];
				
	DtAQSK_out=new int [dim_copiainterna];
	indiceDtAQSK=-1;
	//END copia interna del db Dt_stat.txt 

	PPP=0, NNN=0; 

	/////////////////////////////////////////
	//Per Dt, A e Q	
	double somma2, somma3, sommaquadro2, sommaquadro3;
	double sommacubo, sommacubo2, sommacubo3;
	double sommaquadra, sommaquadra2, sommaquadra3;
	double media2, media3, var2, var3;
	double Skw, Skw2, Skw3, Skw_m3, Skw2_m3, Skw3_m3; // Skewness
	double Kut, Kut2, Kut3, Kut_m4, Kut2_m4, Kut3_m4; // Kurtosis
	
	for(j=0;j<volte;j++){	
		posizione=(int)NumeroRand(dimlog_ripulito);		

		somma=0.0; somma2=0.0; somma3=0.0;
		sommaquadro=0.0; sommaquadro2=0.0; sommaquadro3=0.0;
		sommacubo=0.0; sommacubo2=0.0; sommacubo3=0.0;
		sommaquadra=0.0; sommaquadra2=0.0; sommaquadra3=0.0; 
		
		k=posizione-1;
		for(i=0;i<NumElementi;i++){
			k++;
			if(k==dimlog_ripulito)
				k=0;
			somma+=Dt[k]; sommaquadro+=pow(Dt[k],2); sommacubo+=pow(Dt[k],3);	 sommaquadra+=pow(Dt[k],4);		
			somma2+=sizeA[k]; sommaquadro2+=pow(sizeA[k],2); sommacubo2+=pow(sizeA[k],3); sommaquadra2+=pow(sizeA[k],4);							
			somma3+=sizeQ[k]; sommaquadro3+=pow(sizeQ[k],2); sommacubo3+=pow(sizeQ[k],3); sommaquadra3+=pow(sizeQ[k],4);		

			// per testare il calcolo di Skewness e Kurtosis con matlab stampo i campioni ed il risultato sotto
			//printf("\n%lf",Dt[k]);
		}
		media=somma/NumElementi; media2=somma2/NumElementi; media3=somma3/NumElementi;

		var=sommaquadro+NumElementi*pow(media,2)-2*media*somma; var/=NumElementi;
		var2=sommaquadro2+NumElementi*pow(media2,2)-2*media2*somma2; var2/=NumElementi;
		var3=sommaquadro3+NumElementi*pow(media3,2)-2*media3*somma3; var3/=NumElementi;

		// OCIO devo mettere nella var^(3/2) var^(1.5) oppure var^(3.0/2.0) perchè se metto ^(3/2) il conto lo fa sbagliato!!!
		Skw_m3=sommacubo-3*media*sommaquadro+3*pow(media,2)*somma-NumElementi*pow(media,3); Skw_m3/=NumElementi; Skw=Skw_m3/pow(var,(3.0/2.0)); 		
		Skw2_m3=sommacubo2-3*media2*sommaquadro2+3*pow(media2,2)*somma2-NumElementi*pow(media2,3); Skw2_m3/=NumElementi; Skw2=Skw2_m3/pow(var2,(3.0/2.0));
		Skw3_m3=sommacubo3-3*media3*sommaquadro3+3*pow(media3,2)*somma3-NumElementi*pow(media3,3); Skw3_m3/=NumElementi; Skw3=Skw3_m3/pow(var3,(3.0/2.0));

		Kut_m4=sommaquadra-4*media*sommacubo+6*pow(media,2)*sommaquadro-4*pow(media,3)*somma+NumElementi*pow(media,4); Kut_m4/=NumElementi; Kut=( Kut_m4/(pow(var,2)) )-3;
		Kut2_m4=sommaquadra2-4*media2*sommacubo2+6*pow(media2,2)*sommaquadro2-4*pow(media2,3)*somma2+NumElementi*pow(media2,4); Kut2_m4/=NumElementi; Kut2=( Kut2_m4/(pow(var2,2)) )-3;
		Kut3_m4=sommaquadra3-4*media3*sommacubo3+6*pow(media3,2)*sommaquadro3-4*pow(media3,3)*somma3+NumElementi*pow(media3,4); Kut3_m4/=NumElementi; Kut3=( Kut_m4/(pow(var3,2)) )-3;

		// per testare il calcolo di Skewness e Kurtosis con matlab
		//printf("\n\n\n"); system("pause"); printf("\n%lf\t%lf\t%lf",var,Skw,Kut); system("pause");

		out=0; NNN+=1;
		Out=fopen("DtAQSK_stat.txt","a");
		//fprintf(Out,"\n%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%d",media,media2,media3,var,var2,var3,Skw,Skw2,Skw3,Kut,Kut2,Kut3,out);
		fprintf(Out,"\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d",media,media2,media3,var,var2,var3,Skw,Skw2,Skw3,Kut,Kut2,Kut3,out);
		//fprintf(Out,"\n%.1e;%.1e;%.1e;%d",Skw,Skw2,Skw3,out);
		fclose(Out);	

		if(j==0) Out=fopen("0.txt","w"); else Out=fopen("0.txt","a");		
		//Out=fopen("0.txt","a");		
		fprintf(Out,"\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",media,media2,media3,var,var2,var3,Skw,Skw2,Skw3,Kut,Kut2,Kut3);
		fclose(Out);		

		indiceDtAQSK++;
		Dt_m[indiceDtAQSK]=media;
		A_m[indiceDtAQSK]=media2;
		Q_m[indiceDtAQSK]=media3;
		Dt_v[indiceDtAQSK]=var;
		A_v[indiceDtAQSK]=var2;
		Q_v[indiceDtAQSK]=var3;
		Dt_s[indiceDtAQSK]=Skw;
		A_s[indiceDtAQSK]=Skw2;
		Q_s[indiceDtAQSK]=Skw3;
		Dt_k[indiceDtAQSK]=Kut;
		A_k[indiceDtAQSK]=Kut2;
		Q_k[indiceDtAQSK]=Kut3;
		DtAQSK_out[indiceDtAQSK]=out;
	}
	//END per Dt, A e Q

	// calcolo minimi e massimi per la parte senza tunnel
	// serve per poi racchiudere i punti senza tunnel in un rettangolo e considerare quelli fuori come tunnel
	/*cerca=TrovaIndiceDelMassimo(Dt_m,volte);
	double Dt_mMax=Dt_m[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_m,volte);
	double Dt_mmin=Dt_m[cerca];

	cerca=TrovaIndiceDelMassimo(A_m,volte);
	double A_mMax=A_m[cerca];
	cerca=TrovaIndiceDelMinimo(A_m,volte);
	double A_mmin=A_m[cerca];

	cerca=TrovaIndiceDelMassimo(Q_m,volte);
	double Q_mMax=Q_m[cerca];
	cerca=TrovaIndiceDelMinimo(Q_m,volte);
	double Q_mmin=Q_m[cerca];

	cerca=TrovaIndiceDelMassimo(Dt_v,volte);
	double Dt_vMax=Dt_v[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_v,volte);
	double Dt_vmin=Dt_v[cerca];

	cerca=TrovaIndiceDelMassimo(A_v,volte);
	double A_vMax=A_v[cerca];
	cerca=TrovaIndiceDelMinimo(A_v,volte);
	double A_vmin=A_v[cerca];

	cerca=TrovaIndiceDelMassimo(Q_v,volte);
	double Q_vMax=Q_v[cerca];
	cerca=TrovaIndiceDelMinimo(Q_v,volte);
	double Q_vmin=Q_v[cerca];

	cerca=TrovaIndiceDelMassimo(Dt_s,volte);
	double Dt_sMax=Dt_s[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_s,volte);
	double Dt_smin=Dt_s[cerca];

	cerca=TrovaIndiceDelMassimo(A_s,volte);
	double A_sMax=A_s[cerca];
	cerca=TrovaIndiceDelMinimo(A_s,volte);
	double A_smin=A_s[cerca];

	cerca=TrovaIndiceDelMassimo(Q_s,volte);
	double Q_sMax=Q_s[cerca];
	cerca=TrovaIndiceDelMinimo(Q_s,volte);
	double Q_smin=Q_s[cerca];

	cerca=TrovaIndiceDelMassimo(Dt_k,volte);
	double Dt_kMax=Dt_k[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_k,volte);
	double Dt_kmin=Dt_k[cerca];

	cerca=TrovaIndiceDelMassimo(A_k,volte);
	double A_kMax=A_k[cerca];
	cerca=TrovaIndiceDelMinimo(A_k,volte);
	double A_kmin=A_k[cerca];

	cerca=TrovaIndiceDelMassimo(Q_k,volte);
	double Q_kMax=Q_k[cerca];
	cerca=TrovaIndiceDelMinimo(Q_k,volte);
	double Q_kmin=Q_k[cerca];*/

	//printf("\n\n%lf %lf %lf %lf",Dt_mMax,Dt_mmin,A_mMax,A_mmin); printf("\n\n"); system("pause");
	//END calcolo minimi e massimi per la parte senza tunnel

	//Per DtMix e Amix	
	for(j=0;j<(int)(volte/div_anomaly);j++){	
		posizione=(int)NumeroRand(dimlogMix);		

		somma=0.0; somma2=0.0; somma3=0.0;
		sommaquadro=0.0; sommaquadro2=0.0; sommaquadro3=0.0;
		sommacubo=0.0; sommacubo2=0.0; sommacubo3=0.0;
		sommaquadra=0.0; sommaquadra2=0.0; sommaquadra3=0.0; 
		
		k=posizione-1;
		for(i=0;i<NumElementi;i++){
			k++;
			if(k==dimlogMix)
				k=0;			
			somma+=DtMix[k]; sommaquadro+=pow(DtMix[k],2); sommacubo+=pow(DtMix[k],3); sommaquadra+=pow(DtMix[k],4);				
			somma2+=AMix[k]; sommaquadro2+=pow(AMix[k],2); sommacubo2+=pow(AMix[k],3); sommaquadra2+=pow(AMix[k],4);							
			somma3+=QMix[k]; sommaquadro3+=pow(QMix[k],2); sommacubo3+=pow(QMix[k],3); sommaquadra3+=pow(QMix[k],4);		
		}
		media=somma/NumElementi; media2=somma2/NumElementi; media3=somma3/NumElementi;

		var=sommaquadro+NumElementi*pow(media,2)-2*media*somma; var/=NumElementi;
		var2=sommaquadro2+NumElementi*pow(media2,2)-2*media2*somma2; var2/=NumElementi;
		var3=sommaquadro3+NumElementi*pow(media3,2)-2*media3*somma3; var3/=NumElementi;

		Skw_m3=sommacubo-3*media*sommaquadro+3*pow(media,2)*somma-NumElementi*pow(media,3); Skw_m3/=NumElementi; Skw=Skw_m3/pow(var,(3.0/2.0));
		Skw2_m3=sommacubo2-3*media2*sommaquadro2+3*pow(media2,2)*somma2-NumElementi*pow(media2,3); Skw2_m3/=NumElementi; Skw2=Skw2_m3/pow(var2,(3.0/2.0));
		Skw3_m3=sommacubo3-3*media3*sommaquadro3+3*pow(media3,2)*somma3-NumElementi*pow(media3,3); Skw3_m3/=NumElementi; Skw3=Skw3_m3/pow(var3,(3.0/2.0));

		Kut_m4=sommaquadra-4*media*sommacubo+6*pow(media,2)*sommaquadro-4*pow(media,3)*somma+NumElementi*pow(media,4); Kut_m4/=NumElementi; Kut=( Kut_m4/(pow(var,2)) )-3;
		Kut2_m4=sommaquadra2-4*media2*sommacubo2+6*pow(media2,2)*sommaquadro2-4*pow(media2,3)*somma2+NumElementi*pow(media2,4); Kut2_m4/=NumElementi; Kut2=( Kut2_m4/(pow(var2,2)) )-3;
		Kut3_m4=sommaquadra3-4*media3*sommacubo3+6*pow(media3,2)*sommaquadro3-4*pow(media3,3)*somma3+NumElementi*pow(media3,4); Kut3_m4/=NumElementi; Kut3=( Kut_m4/(pow(var3,2)) )-3;

		out=1; PPP+=1;
		Out=fopen("DtAQSK_stat.txt","a");		
		//fprintf(Out,"\n%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%d",media,media2,media3,var,var2,var3,Skw,Skw2,Skw3,Kut,Kut2,Kut3,out);
		fprintf(Out,"\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d",media,media2,media3,var,var2,var3,Skw,Skw2,Skw3,Kut,Kut2,Kut3,out);
		//fprintf(Out,"\n%.1e;%.1e;%.1e;%d",Skw,Skw2,Skw3,out);
		fclose(Out);				

		if(j==0) Out=fopen("1.txt","w"); else Out=fopen("1.txt","a");		
		//Out=fopen("1.txt","a");
		fprintf(Out,"\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",media,media2,media3,var,var2,var3,Skw,Skw2,Skw3,Kut,Kut2,Kut3);
		fclose(Out);		

		indiceDtAQSK++;
		Dt_m[indiceDtAQSK]=media;
		A_m[indiceDtAQSK]=media2;
		Q_m[indiceDtAQSK]=media3;
		Dt_v[indiceDtAQSK]=var;
		A_v[indiceDtAQSK]=var2;
		Q_v[indiceDtAQSK]=var3;
		Dt_s[indiceDtAQSK]=Skw;
		A_s[indiceDtAQSK]=Skw2;
		Q_s[indiceDtAQSK]=Skw3;
		Dt_k[indiceDtAQSK]=Kut;
		A_k[indiceDtAQSK]=Kut2;
		Q_k[indiceDtAQSK]=Kut3;
		DtAQSK_out[indiceDtAQSK]=out;
	}
	//END per DtMix e Amix
	/////////////////////////////////////////
	//return;

	//END per db e dbMix piglio a caso un punto di partenza e piglio da lì 5000 campioni e ne calcolo media e varianza, tutto ciò ripetutamente
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------

	// calcolo minimi e massimi per tutto il db con e senza tunnel
	/*volte=indiceDtAQSK+1;
	cerca=TrovaIndiceDelMassimo(Dt_m,volte);
	Dt_mMax2=Dt_m[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_m,volte);
	Dt_mmin2=Dt_m[cerca];

	cerca=TrovaIndiceDelMassimo(A_m,volte);
	A_mMax2=A_m[cerca];
	cerca=TrovaIndiceDelMinimo(A_m,volte);
	A_mmin2=A_m[cerca];

	cerca=TrovaIndiceDelMassimo(Q_m,volte);
	Q_mMax2=Q_m[cerca];
	cerca=TrovaIndiceDelMinimo(Q_m,volte);
	Q_mmin2=Q_m[cerca];

	cerca=TrovaIndiceDelMassimo(Dt_v,volte);
	Dt_vMax2=Dt_v[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_v,volte);
	Dt_vmin2=Dt_v[cerca];

	cerca=TrovaIndiceDelMassimo(A_v,volte);
	A_vMax2=A_v[cerca];
	cerca=TrovaIndiceDelMinimo(A_v,volte);
	A_vmin2=A_v[cerca];

	cerca=TrovaIndiceDelMassimo(Q_v,volte);
	Q_vMax2=Q_v[cerca];
	cerca=TrovaIndiceDelMinimo(Q_v,volte);
	Q_vmin2=Q_v[cerca];

	cerca=TrovaIndiceDelMassimo(Dt_s,volte);
	Dt_sMax2=Dt_s[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_s,volte);
	Dt_smin2=Dt_s[cerca];

	cerca=TrovaIndiceDelMassimo(A_s,volte);
	A_sMax2=A_s[cerca];
	cerca=TrovaIndiceDelMinimo(A_s,volte);
	A_smin2=A_s[cerca];

	cerca=TrovaIndiceDelMassimo(Q_s,volte);
	Q_sMax2=Q_s[cerca];
	cerca=TrovaIndiceDelMinimo(Q_s,volte);
	Q_smin2=Q_s[cerca];

	cerca=TrovaIndiceDelMassimo(Dt_k,volte);
	Dt_kMax2=Dt_k[cerca];
	cerca=TrovaIndiceDelMinimo(Dt_k,volte);
	Dt_kmin2=Dt_k[cerca];

	cerca=TrovaIndiceDelMassimo(A_k,volte);
	A_kMax2=A_k[cerca];
	cerca=TrovaIndiceDelMinimo(A_k,volte);
	A_kmin2=A_k[cerca];

	cerca=TrovaIndiceDelMassimo(Q_k,volte);
	Q_kMax2=Q_k[cerca];
	cerca=TrovaIndiceDelMinimo(Q_k,volte);
	Q_kmin2=Q_k[cerca];*/
	// END calcolo minimi e massimi per tutto il db con e senza tunnel

	//return;

	printf("\n\n---------------------------------------------------");
	printf("\n\tErrori classificatori\n"); 
	// classificatore Lineare (dal matlab)
	double errLineare=0;
	int out_classificatoreLineare;
	//printf("\n\n");
	for(j=0;j<=0;j++){ // j<=indiceDtAQSK

		out_classificatoreLineare=ClassificatoreLineare(Dt_m[j], A_m[j]);
		if(DtAQSK_out[j]!=out_classificatoreLineare)
			errLineare++;
		
		//printf("%f;%f;%d\t%d\t%d\n",Dt_m[j],A_m[j],DtAQSK_out[j],out_classificatoreLineare);		
	}
	errLineare/=indiceDtAQSK;	
	//printf("\nLineare=%.1f%%",errLineare*100); 
	//END Lineare

	// Lineare3D
	errLineare=0;
	//printf("\n\n");
	for(j=0;j<=0;j++){ // j<=indiceDtAQSK

		out_classificatoreLineare=ClassificatoreLineare3D(Dt_m[j], A_m[j], Q_m[j]);		
		if(DtAQSK_out[j]!=out_classificatoreLineare)
			errLineare++;					
	}
	errLineare/=indiceDtAQSK	;
	//printf("\n\tLineare3D=%.0f%%",errLineare*100); 	
	//printf("\n\tLineare3D=%lf",errLineare); 	
	//END Lineare3D	

	// Lineare 6D
	errLineare=0;
	//printf("\n\n");
	for(j=0;j<=-1;j++){ //j<=indiceDtAQSK;
		
		out_classificatoreLineare=ClassificatoreLineare6D(Dt_m[j], A_m[j], Q_m[j], Dt_v[j], A_v[j], Q_v[j]);
		if(DtAQSK_out[j]!=out_classificatoreLineare)
			errLineare++;					
	}
	errLineare/=indiceDtAQSK	;	
	//printf("\n\tLineare6D=%.0f%%\t%e",errLineare*100,errLineare); 	
	//END Lineare 6D	

	

	

	double pesi[4]; 
	for(j=0;j<4;j++) 
		pesi[j]=0.001;	
	//pesi[0]=0; pesi[1]=0; pesi[2]=0; pesi[3]=1; 
	//Voting12DCFSQP(pesi);

	/*printf("\n\nPesi dopo l'ottimizzazione");
	for(j=0;j<4;j++)
		printf("\npesi[%d]=%f",j,pesi[j]);
	system("pause");*/

	/*for(pesi[0]=0.01;pesi[0]<=0.9;pesi[0]+=0.025)
		for(pesi[1]=0.01;pesi[1]<=0.9;pesi[1]+=0.025)
			for(pesi[2]=0.01;pesi[2]<=0.9;pesi[2]+=0.025)
				for(pesi[3]=0.01;pesi[3]<=0.9;pesi[3]+=0.025){

					printf("\n\nPesi estratti");	
					for(j=0;j<4;j++)
						printf("\npesi[%d]=%f",j,pesi[j]);
					
					errore=ErroreVoting12DD(indiceDtAQSK,pesi);
					printf("\n\tVoting12D=%.0f%%\t%e",(errore)*100,errore);
					if(errore<0.35)
						system("pause");
				}*/
		
	//system("pause");

	//printf("\n\tVoting12D=%.0f%%\t%e",(ErroreVoting12DD(indiceDtAQSK,pesi))*100,(ErroreVoting12DD(indiceDtAQSK,pesi))); 	

	// Rettangolo
	/*double errRettangolo=0;
	int out_classificatoreRettangolo; 
	//printf("\n\n");
	for(j=0;j<=0;j++){ // j<=indiceDtAQSK

		out_classificatoreRettangolo=ClassificatoreRettangolo(Dt_mmin, A_mmin, Dt_mMax, A_mMax, Dt_m[j], A_m[j]);
		if(DtAQSK_out[j]!=out_classificatoreRettangolo)
			errRettangolo++;		
		
		//printf("%f;%f;%d\t%d\n",Dt_m[j],A_m[j],DtAQSK_out[j],out_classificatoreRettangolo);		
	}
	errRettangolo/=indiceDtAQSK;		
	//printf("\nRettangolo=%.0f%%",errRettangolo*100); 
	//printf("\nLimiti rettangolo: %lf %lf %lf %lf\n\n\t",Dt_mMax,Dt_mmin,A_mMax,A_mmin);
	*/
	//END Rettangolo

	

	//printf("\nLimiti rettangolo: %lf %lf %lf %lf\n\n\t",Dt_mMax,Dt_mmin,A_mMax,A_mmin);

	////////////////////////////////////////// stampo tutto il dataset //////////////////////////////////////////
	/*Out=fopen("Dataset.txt","w");
	fprintf(Out,"\nQ\tA\tDt");
	for(i=0;i<dimlog;i++){
		if(sizeA[i]>0.0)
			fprintf(Out,"\n%lf\t%lf\t%lf",sizeQ[i],sizeA[i],tA[i]-tQ[i]);
	}
	fclose(Out);
	return;*/
	////////////////////////////////////////// END stampo tutto il dataset //////////////////////////////////////////

	////////////////////////////////////////// stampo il dataset di un IP specifico //////////////////////////////////////////		
	/*Out=fopen("150.145.0.118.txt","w");
	fprintf(Out,"\nQ\tA\tDt");
	char *stringadacercare2="150.145.0.118";	
	for(i=0;i<dimlog;i++){
		if ( strncmp(src[i],stringadacercare2,(unsigned)strlen(stringadacercare2)) == 0){
			  //printf ("\n%s al posto %d",src[i],i);

			  fprintf(Out,"\n%lf\t%lf\t%lf",sizeQ[i],sizeA[i],tA[i]-tQ[i]);
		}
	}	
	fclose(Out);
	return;*/
	////////////////////////////////////////// END stampo il dataset di un IP specifico //////////////////////////////////////////

	////////////////////////////////////////// Stampo tanti file per ogni IP src //////////////////////////////////////////	
	int num=-1;
	int indice=0;
	char **estratti;
	estratti=new char *[dimlog];
	for(p=0;p<dimlog;p++)
		estratti[p]=new char[20];	
	/*int Elementiin_estratti;	
	bool giapresente;		
	
	for(i=0;i<dimlog;i++){

		giapresente=false;
		
		strncpy (local,src[i],20);
		for(j=0;j<num;j++){
			if ( strncmp(estratti[j],local,(unsigned)strlen(local)) == 0)
				giapresente=true;
		}

		if(!giapresente){

			num++;
			strncpy (estratti[num],local,20);
			
			Elementiin_estratti=0;
			for(j=0;j<dimlog;j++){
				if ( strncmp(src[j],local,(unsigned)strlen(local)) == 0)
					 Elementiin_estratti++;
			}			

			if(Elementiin_estratti>10000){			
				indice++;
				Out=fopen(local,"w");
				fprintf(Out,"\nQ%d\tA%d\tDt%d",indice,indice,indice);
				for(j=0;j<dimlog;j++){
					if ( (strncmp(src[j],local,(unsigned)strlen(local)) == 0) && (sizeA[j]>0.0) )
						 fprintf(Out,"\n%lf\t%lf\t%lf",sizeQ[j],sizeA[j],tA[j]-tQ[j]);
				}
				fclose(Out);
			}
		
		}//END if(!giapresente)
	
	}//END for(i=0;i<dimlog;i++)*/
	////////////////////////////////////////// END Stampo tanti file per ogni IP src //////////////////////////////////////////

	fclose(dns);	
	fclose(dns2);	

	delete sizeQ;
	delete sizeA;
	delete tQ;
	delete tA;
	delete N_filelog;	
	delete IdDns;
	delete Dt;
	for(i=0;i<dimlog;i++)
			delete( estratti[i] );
	delete estratti;
	for(i=0;i<dimlog;i++)
			delete( src[i] );
	delete src;
	for(i=0;i<dimlog;i++)
			delete( dst[i] );
	delete dst;

	delete sizeQ2;
	delete sizeA2;
	delete tQ2;
	delete tA2;
	delete N_filelog2;	
	delete IdDns2;
	delete Dt2;
	for(i=0;i<dimlog2;i++)
			delete( src2[i] );
	delete src2;
	for(i=0;i<dimlog2;i++)
			delete( dst2[i] );
	delete dst2;

	delete DtMix;
	delete AMix;
	delete QMix;

	delete Dt_stat_m;
	delete Dt_stat_v;
	delete Dt_stat_out;

	delete Dt_m;
	delete A_m;
	delete Q_m;
	delete Dt_v;
	delete A_v;
	delete Q_v;
	delete Dt_s;
	delete A_s;
	delete Q_s;
	delete Dt_k;
	delete A_k;
	delete Q_k;
	delete DtAQSK_out;	

	printf("\n\n"); system("pause");
	
}//fine main
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void Permuta(int volte, int Dim, double *x, double *y, double *z){

	double local;
	int posizione1, posizione2;
	
	for(int i=0; i<volte; i++){

		posizione1=(int)NumeroRand(Dim); if(posizione1>=Dim) posizione1=(int)NumeroRand(Dim);
		if(posizione1>=Dim){printf("\nposizione1>=Dim"); system("pause");}
	   
		posizione2=(int)NumeroRand(Dim); if(posizione2>=Dim) posizione2=(int)NumeroRand(Dim);
		if(posizione2>=Dim){printf("\nposizione2>=Dim"); system("pause");}
	
		local=x[posizione1];
		x[posizione1]=x[posizione2];
		x[posizione2]=local;

		local=y[posizione1];
		y[posizione1]=y[posizione2];
		y[posizione2]=local;

		local=z[posizione1];
		z[posizione1]=z[posizione2];
		z[posizione2]=local;

	}

}
////////////////////////////////////////////////////////////////////

int ClassificatoreLineare(double x, double y){

	//coefficienti del matlab
	double K= 20.6363;
	double L1=  -6.5368;
	double L2= -0.0728;

	if( (K+L1*x+L2*y) <= 0 )
		return 1;
	else
		return 0;
}

int ClassificatoreLineare3D(double x, double y, double z){
	
	/*double K= 240.4851; // ssh ieiit 5e3 1e4 10% => err=10%
	double L1= 1.5440;
	double L2= -0.3376;
	double L3= -1.7219;*/
	/*double K= 301.0462; // come su, ma con area2 => err=8%
	double L1= -17.0014;
	double L2= -0.3623;
	double L3= -2.4835;*/
	/*double K= 50.0556; // come su, ma con ieiit+area2 => err=15%
	double L1= 6.5319;
	double L2= -0.1403;
	double L3= -0.1838;*/

	double K= 50.0556; // come su, ma con ieiit+area2 => err=15%
	double L1= 6.5319;
	double L2= -0.1403;
	double L3= -0.1838;

	if( (K+L1*x+L2*y+L3*z) <= 0 )
		return 1;
	else
		return 0;
}

int ClassificatoreLineare6D(double x, double y, double z, double x2, double y2, double z2){
	
	/* coefficienti del matlab per il ssh-ieiit, solo medie, mix 1%, 10^4	*/
	double K= 56.0153;
	double L1= 1.6974;
	double L2= -0.2210;
	double L3= -0.1629; 
	double L4= -1.6686e-004;
	double L5= 4.4283e-004;
	double L6= 3.2641e-004; 

	if( (K+L1*x+L2*y+L3*z+L4*x2+L5*y2+L6*z2) <= 0 )
		return 1;
	else
		return 0;
}

////////////////////////////////////////////////////////////////////
double ErroreLineare12D_ripetuto_amaggioranza(int Size, int classificatore, int num_campioni){ 

	double errLineare=0;
	int i=0, j, out_classificatoreLineare, conta0, conta1;
	//printf("\n\n");	

	while(i<=Size){ // la prima metà sono no tunnel output previsto 0

		conta0=0; conta1=0;
		for(j=i;j<i+num_campioni;j++){
			out_classificatoreLineare=ClassificatoreLineare12D(Dt_m[j], A_m[j], Q_m[j], Dt_v[j], A_v[j], Q_v[j], Dt_s[j], A_s[j], Q_s[j], Dt_k[j], A_k[j], Q_k[j],classificatore);
			if(out_classificatoreLineare==0)
				conta0++;
			if(out_classificatoreLineare==1)
				conta1++;
		}
		if(conta0>=conta1)
			out_classificatoreLineare=0;
		else
			out_classificatoreLineare=1;
				
		if(i<Size/2){
			if(0!=out_classificatoreLineare)
				errLineare++;		
		}else{
			if(1!=out_classificatoreLineare)
				errLineare++;	
		}

		i=j;
	}

	errLineare/=Size;	
	
	return errLineare;
}
////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////
double ErroreLineare12D(int Size, int classificatore, double *perf){ 

	double errLineare=0, FP=0, TP=0, TN=0, FN=0;
	int j, out_classificatoreLineare;
	//printf("\n\n");
	for(j=0;j<=Size;j++){ 
		
		out_classificatoreLineare=ClassificatoreLineare12D(Dt_m[j], A_m[j], Q_m[j], Dt_v[j], A_v[j], Q_v[j], Dt_s[j], A_s[j], Q_s[j], Dt_k[j], A_k[j], Q_k[j],classificatore);
		if(DtAQSK_out[j]!=out_classificatoreLineare)
			errLineare++;					

		if(j<Size/2){ // la prima metà sono campioni clean (ho preso questo if da 'ErroreNeurale12D_ripetuto_amaggioranza()')
			if(0!=out_classificatoreLineare) {FP++;}	else TN++; 							
		}else{ // la seconda metà sono campioni di tunnel
			if(1!=out_classificatoreLineare) {FN++;} else TP++; 					
		}

	}
	errLineare/=Size;	
	FP/=Size; TP/=Size; TN/=Size; FN/=Size;

	if(perf) { perf[0]=errLineare; perf[1]=1-errLineare; perf[2]=(TP)/(TP+FP); perf[3]=FP; perf[4]=FN; perf[5]=(TP+TN)/(TP+TN+FP+FN); }
	
	return errLineare;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double ErroreVoting12DD(int Size, double *pesi){

	double errLineare=0;
	int j, out_classificatoreLineare;
	//printf("\n\n");
	for(j=0;j<=Size;j++){ 
		
		out_classificatoreLineare=VotingLineare12D(Dt_m[j], A_m[j], Q_m[j], Dt_v[j], A_v[j], Q_v[j], Dt_s[j], A_s[j], Q_s[j], Dt_k[j], A_k[j], Q_k[j], pesi);
		if(DtAQSK_out[j]!=out_classificatoreLineare)
			errLineare++;					
	}
	errLineare/=Size;	
	
	return errLineare;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int VotingLineare12D(double x, double y, double z, double x2, double y2, double z2, 
					 double x3, double y3, double z3, double x4, double y4, double z4, double *pesi){

	 double somma1=0, somma0=0;
	 int i;

	 for(i=1;i<=4;i++){

		 if(ClassificatoreLineare12D(x,y,z,x2,y2,z2,x3,y3,z3,x4,y4,z4, i)==1)
			 somma1+=pesi[i];
		 else
			 somma0+=pesi[i];

	 }//END for(i=1;i<=4;i++)

	 if(somma1>=somma0)
		 return 1;
	 else
		return 0;

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double mapminmax(double x, double xmin, double xMax){

	return ( (ymax-ymin)*(x-xmin)/(xMax-xmin) + ymin );
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double mapminmax_reverse(double y, double xmin, double xMax){
	
	return ( (y-ymin)*(xMax-xmin)/(ymax-ymin)+xmin );
}
////////////////////////////////////////////////////////////////////



int ClassificatoreLineare12D(double x, double y, double z, double x2, double y2, double z2, 
												double x3, double y3, double z3, double x4, double y4, double z4, int classificatore){
			
	double K, L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12;	

	/* coefficienti del matlab per il ieiit, m/v/s/k, mix 1%, 10^4, srand=5
	classificatore=1 --> p2p (seguendo la colonna del 3.tentativo.doc)
	classificatore=2 --> ssh
	classificatore=3 --> wget
	classificatore=4 --> mix (all)
	*/
	if(classificatore==1){
		K= 1.2170e+003;
		L1= 12.0812;
		L2= -0.8598;
		L3= -11.7933; 
		L4= -0.0011;
		L5= 0.0019;
		L6= 0.0096; 
		L7= -0.1075;
		L8= -8.3843;
		L9= 0.0448; 
		L10= 7.5261e-004;
		L11= -3.1261;
		L12= -3.4372e-007; 
	}
	if(classificatore==2){
		
		if(0){
			K= 27.7303;	// ieiit	
  L1= -0.1667;
  L2= -0.0280;
   L3=-0.5580;
    L4=0.0000;
    L5=0.0002;
   L6= 0.0005;
  L7= -0.3671;
   L8=19.3500;
   L9=-0.0159;
    L10=0.0099;
  L11= -2.4458;
  L12= -0.0000;
		}

		if(1){
			K= 145.1108; // area
		 L1=-23.0898;
    L2=-0.0699;
  L3= -1.5359;
   L4= 0.7961;
  L5= -0.0003;
    L6=0.0032;
  L7=  -0.0371;
   L8= 7.6654;
  L9= -0.0313;
     L10=0.0007;
 L11=  -0.4904;
   L12= 1.2877;
		}

	}	
	if(classificatore==3){
		K= 12.8757;
		L1= 1.0920;
		L2= -0.0379; 
		L3= -0.5727; 
		L4= -9.3261e-005;
		L5= 2.0279e-004;
		L6= 7.2557e-004; 
		L7= -0.3757;
		L8= 28.6215;
		L9= -0.0033; 
		L10= 0.0036;
		L11= -3.3188;
		L12= -1.2266e-007; 
	}	
	if(classificatore==4){
		K= 210.2321;

		/*L1=1.6066; // ieiit
   L2=-0.2114;
   L3=-3.3545;
   L4=-0.0001;
   L5= 0.0006;
    L6=0.0025;
   L7=-0.4460;
   L8=11.6225;
   L9=-0.0179;
    L10=0.0093;
   L11=-1.8507;
    L12=0.0000;*/

	 L1=2.9357; // area
    L2=-0.1023;
   L3=-2.1850;
   L4=-0.1725;
   L5=-0.0002;
    L6=0.0065;
  L7= -0.0432;
   L8= 8.4426;
   L9= -0.6542;
    L10=0.0008;
  L11= -0.6379;
   L12= 1.8991;
	}	

	if( (K+L1*x+L2*y+L3*z+L4*x2+L5*y2+L6*z2	+L7*x3+L8*y3+L9*z3+L10*x4+L11*y4+L12*z4) <= 0 )
		return 1;
	else
		return 0;
}


int ClassificatoreQ(double m, double v){

	if( (-0.0696723+465.548*m+-18.5204*v+-12529.8*m*m+1327.17*m*v+-38.4875*v*v)<= 0)
		return 1;
	else
	return 0;
}

int ClassificatoreRettangolo(double Dt_mmin, double A_mmin, double Dt_mMax, double A_mMax, double Dt_mSample,  double A_mSample){

	int out=0;

	if(Dt_mSample<Dt_mmin)
		out=1;
	if(Dt_mSample>Dt_mMax)
		out=1;

	if(A_mSample<A_mmin)
		out=1;
	if(A_mSample>A_mMax)
		out=1;

	return out;
}



int ClassificatoreRettangolo3D(double Dt_mmin, double A_mmin, double Q_mmin, double Dt_mMax, double A_mMax, double Q_mMax, 
							   double Dt_mSample,  double A_mSample,  double Q_mSample){

	int out=0;

	if(Dt_mSample<Dt_mmin)
		out=1;
	if(Dt_mSample>Dt_mMax)
		out=1;

	if(A_mSample<A_mmin)
		out=1;
	if(A_mSample>A_mMax)
		out=1;

	if(Q_mSample<Q_mmin)
		out=1;
	if(Q_mSample>Q_mMax)
		out=1;

	return out;
}

int ClassificatoreRettangolo6D(double Dt_mmin, double A_mmin, double Q_mmin, double Dt_mMax, double A_mMax, double Q_mMax,
							   double Dt_vmin, double A_vmin, double Q_vmin, double Dt_vMax, double A_vMax, double Q_vMax,
							   double Dt_mSample,  double A_mSample,  double Q_mSample,
							   double Dt_vSample,  double A_vSample,  double Q_vSample){

	int out=0;

	if(Dt_mSample<Dt_mmin)
		out=1;
	if(Dt_mSample>Dt_mMax)
		out=1;
	if(Dt_vSample<Dt_vmin)
		out=1;
	if(Dt_vSample>Dt_vMax)
		out=1;

	if(A_mSample<A_mmin)
		out=1;
	if(A_mSample>A_mMax)
		out=1;
	if(A_vSample<A_vmin)
		out=1;
	if(A_vSample>A_vMax)
		out=1;

	if(Q_mSample<Q_mmin)
		out=1;
	if(Q_mSample>Q_mMax)
		out=1;
	if(Q_vSample<Q_vmin)
		out=1;
	if(Q_vSample>Q_vMax)
		out=1;

	return out;
}

int ClassificatoreRettangolo12D(double Dt_mmin, double A_mmin, double Q_mmin, double Dt_mMax, double A_mMax, double Q_mMax,
							   double Dt_vmin, double A_vmin, double Q_vmin, double Dt_vMax, double A_vMax, double Q_vMax,
							   double Dt_smin, double A_smin, double Q_smin, double Dt_sMax, double A_sMax, double Q_sMax,
							   double Dt_kmin, double A_kmin, double Q_kmin, double Dt_kMax, double A_kMax, double Q_kMax,
							   double Dt_mSample,  double A_mSample,  double Q_mSample,
							   double Dt_vSample,  double A_vSample,  double Q_vSample,
							   double Dt_sSample,  double A_sSample,  double Q_sSample,
							   double Dt_kSample,  double A_kSample,  double Q_kSample){

	int out=0;

	if(Dt_mSample<Dt_mmin)
		out=1;
	if(Dt_mSample>Dt_mMax)
		out=1;
	if(Dt_vSample<Dt_vmin)
		out=1;
	if(Dt_vSample>Dt_vMax)
		out=1;
	if(Dt_sSample<Dt_smin)
		out=1;
	if(Dt_sSample>Dt_sMax)
		out=1;
	if(Dt_kSample<Dt_kmin)
		out=1;
	if(Dt_kSample>Dt_kMax)
		out=1;

	if(A_mSample<A_mmin)
		out=1;
	if(A_mSample>A_mMax)
		out=1;
	if(A_vSample<A_vmin)
		out=1;
	if(A_vSample>A_vMax)
		out=1;
	if(A_sSample<A_smin)
		out=1;
	if(A_sSample>A_sMax)
		out=1;
	if(A_kSample<A_kmin)
		out=1;
	if(A_kSample>A_kMax)
		out=1;

	if(Q_mSample<Q_mmin)
		out=1;
	if(Q_mSample>Q_mMax)
		out=1;
	if(Q_vSample<Q_vmin)
		out=1;
	if(Q_vSample>Q_vMax)
		out=1;
	if(Q_sSample<Q_smin)
		out=1;
	if(Q_sSample>Q_sMax)
		out=1;
	if(Q_kSample<Q_kmin)
		out=1;
	if(Q_kSample>Q_kMax)
		out=1;

	return out;
}


////////////////////////////////////////////////////////////////////
double NumeroRand(double max){

	double out;
	
	double NumeroRand0ToRandMax=rand();
	out=NumeroRand0ToRandMax/RAND_MAX;
	out=out*max;

	return out;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double Segno(double numero){

	if(numero<0)
		return (-1);
	else
		return (+1);
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int TrovaIndiceDelMinimo(double *Array, int MaxDim){

double minimo=1000000000000.0;
int indice;
for(int i=0;i<MaxDim;i++)
	if(Array[i]<minimo){
                minimo=Array[i];
				indice=i;
	}
				
return indice;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double TrovaMinimo(double *Array, int MaxDim){

double minimo=1000000000000.0;
for(int i=0;i<MaxDim;i++)
        if(Array[i]<minimo)
                minimo=Array[i];
return minimo;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double TrovaMassimo(double *Array, int MaxDim){

double Max = -1000000000000.0;
for(int i=0;i<MaxDim;i++)
        if(Array[i]>Max)
                Max=Array[i];
return Max;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int TrovaIndiceDelMassimo(double *Array, int MaxDim){

double Max = -1000000000000.0;
int indice;
for(int i=0;i<MaxDim;i++)
	if(Array[i]>Max){
                Max=Array[i];
				indice=i;
	}
				
return indice;
}
////////////////////////////////////////////////////////////////////