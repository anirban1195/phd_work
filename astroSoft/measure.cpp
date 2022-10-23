#include <fitsio.h>
#include <fitsio2.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979328462643

int main(int argc, char **argv) {

    fitsfile *faptr;
    long naxes[2];
    int nfound;
    int anynull;
    float nullval;
    char keyname[4096];
    char comment[4096];
    char value[4096];
    int status;
    float *array;
    double *data;
    float *data2;

    status=0;
    if (fits_open_file(&faptr,argv[1],READONLY,&status)) {
        printf("Error:  Could not open %s\n",argv[1]);
        exit(1);
    }
    fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
    array=static_cast<float*>(calloc(naxes[0]*naxes[1],sizeof(float)));
    fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,array,&anynull,&status);
    fits_close_file(faptr,&status);

    double tot2 = 0.0;
    for (long i=0;i<naxes[0];i++) {
        for (long j=0;j<naxes[1];j++) {
            tot2+=array[i*naxes[1]+j];
        }
    }
    //printf("Total counts in the image:  %lf\n",tot2);

    int cut=50;
    /*if(naxes[0] <= naxes[1])
    {
	cut = naxes[0]/2;
    }
    else
    {
	cut = naxes[1]/2;
    }

    if (cut < 25)
    {
	cut = 25;
    }
    if (cut > 50)
    {
	cut = 50;
    }*/

    printf("Cut:  %d \n",cut);

    if ((naxes[0] > (2*cut)) || (naxes[1] > (2*cut))) {
            int maxx; int maxy; float maxv = 0.0;
               for (int i = 0; i < naxes[0]; i++) {
                   for (int j=0; j < naxes[1]; j++) {
                       if (array[j*naxes[0]+i] > maxv) {
                           maxv=array[j*naxes[0]+i];
                           maxx=i;
                           maxy=j;
                       }
                   }
               }
	       printf("Maxx, maxy, maxv:  %d   %d %f \n",maxx, maxy, maxv);
               int lowx=maxx-cut;
               if (lowx < 0) lowx=0;
               int lowy=maxy-cut;
               if (lowy < 0) lowy=0;
               int highx = maxx+cut;
               if (highx > naxes[0]-1) highx=naxes[0]-1;
               int highy = maxy+cut;
               if (highy > naxes[1]-1) highy=naxes[1]-1;
		printf("X ,y :  %d   %d %d %d \n",highx, lowx, highy, lowy);
               //printf("Image too large: using a cut out of (%d:%d,%d:%d)\n",lowx,highx,lowy,highy);
               data2=static_cast<float*>(calloc((highx-lowx)*(highy-lowy),sizeof(float)));
               for (int i = lowx; i < highx; i++) {
                   for (int j = lowy; j< highy; j++) {
                       data2[(j-lowy)*(highx-lowx)+(i-lowx)]=array[j*naxes[0]+i];
                   }
               }
               naxes[0]=highx-lowx;
               naxes[1]=highy-lowy;
               for (int i = 0; i < naxes[0]; i++) {
                   for (int j=0; j < naxes[1]; j++) {
                       array[j*naxes[0]+i]=data2[j*naxes[0]+i];
                   }
               }
    }

    double tot3 = 0.0;
    for (long i=0;i<naxes[0];i++) {
        for (long j=0;j<naxes[1];j++) {
            tot3+=array[j*naxes[0]+i];
        }
    }
    printf("Total counts in the image:  %lf\n",tot3);

    double sizex=static_cast<double>(naxes[0]);
    double sizey=static_cast<double>(naxes[1]);
    data=static_cast<double*>(calloc(naxes[0]*naxes[1],sizeof(double)));
    for (int i = 0; i < naxes[0]; i++) {
        for (int j=0; j < naxes[1]; j++) {
            data[j*naxes[0]+i]=static_cast<double>(array[j*naxes[0]+i]);
        }
    }

    long nn=naxes[0]*naxes[1];
    for(int i=0;i<nn-1;i++) {
        for(int j=0;j<nn-i-1;j++) {
            if(array[j]>array[j+1]) {
                double t;
                t=array[j];
                array[j]=array[j+1];
                array[j+1]=t;
            }
        }
    }
    long n0 = (nn+1) / 2 - 1;
    long med = array[n0];
    //printf("Median of image is %ld counts\n", med);
    //printf("\n");
    long count=0;
    if (med == 0) {
        for (int i=0;i<sizex;i++) {
            for (int j=0;j<sizey;j++) {
                if (data[j*naxes[0]+i] != 0) {
                    count++;
                }
            }
        }
    } else {
        count=naxes[0]*naxes[1];
    }
    double *x, *y, *w, *weight;
    long xx, yy;
    x=static_cast<double*>(calloc(count,sizeof(double)));
    y=static_cast<double*>(calloc(count,sizeof(double)));
    w=static_cast<double*>(calloc(count,sizeof(double)));
    weight=static_cast<double*>(calloc(count,sizeof(double)));
    count=0;
    for (int i=0;i<naxes[0];i++) {
        for (int j=0;j<naxes[1];j++) {
            if (data[j*naxes[0]+i] != 0) {
            x[count]=j-(sizey/2.0-0.5);
            y[count]=i-(sizex/2.0-0.5);
            w[count]=data[j*naxes[0]+i];
            count++;
            }
        }
    }

    double total = 0.0;
    for (long i=0;i<count;i++) total+=w[i];
    FILE *fp;
    fp = fopen("/media/anirban/anirban/cpp_temp/test.txt", "w+");
    double a=0.0;
    double b=0.0;
    double mux=0.0;
    double muy=0.0;
    double alphax=6;
    double alphax1 = 6;
    double alphay=6;
    double alphay1 = 6;
    double alphaxy=0.0;
    double background=med;
    double t1, t2, t3, t4, t5, t6, t7;
    double sigmaxx, sigmayy, sigmaxy, rms, signal, e1, e2, prevsigmaxx, prevsigmayy, prevsigmaxy;
    int converged = 0;
    double ellip, pa;
    ellip=pa=0.0;
    for (long iteration = 0; iteration < 100; iteration++) {
	//Catching bad things before happening
	if(isnan(sigmaxx) == 1 or isnan(sigmayy) == 1 or isnan(sigmaxy) == 1 or isnan(signal) == 1 or isinf(sigmaxx) == 1 or isinf(sigmayy) == 1 or isinf(sigmaxy) == 1 or isinf(signal) == 1)
	{
		converged =0; 
		break;
	}
	if(abs(sigmaxx) >1000  or abs(sigmayy) > 1000 or abs(sigmaxy) > 1000)
	{
		converged =0; 
		break;
	}
        t1=0.0;
        t2=0.0;
        t3=0.0;
        t4=0.0;
        t5=0.0;
        t6=0.0;
        t7=0.0;
	alphax1 = alphax;
	alphay1 = alphay;
        //If gaussian is too sharply peaked restrict weight alphas
	if(alphax1 < 0.75)
	    alphax1 = 0.75;
	if(alphay1 < 0.75)
	    alphay1 = 0.75;
	// Making sure the sqrt in weight isnt giving nans if root of negative number
	while(pow((alphaxy/alphax1/alphay1),2.0) >1.0 )
	{
		if(alphaxy > 0)
		alphaxy = alphaxy - 0.1012345;
		if(alphaxy < 0)
		alphaxy = alphaxy + 0.1012345;
	}
        for (long i = 0; i < count; i++) {
            weight[i]=(exp(-(pow((x[i]-mux),2.0)/alphax1/alphax1-
                             2.0*alphaxy/alphax1/alphax1/alphay1/alphay1*(x[i]-mux)*(y[i]-muy)+
                             pow((y[i]-muy),2.0)/alphay1/alphay1)/2.0/(1-pow((alphaxy/alphax1/alphay1),2.0))))/
                (2*PI*alphax1*alphay1*sqrt(1-pow((alphaxy/alphax1/alphay1),2.0)));
	    
	    //Testing stuff
	    //if(isnan(weight[i]))
	    //{
		//a=  sqrt(1-pow((alphaxy/alphax1/alphay1),2.0));
		//b = pow((alphaxy/alphax1/alphay1),2.0);
		//printf("%lf %lf %lf %lf %lf %lf \n", weight[i], a, b, alphaxy, alphax1, alphay1);
	    //}
	    //if(weight[i] < 0.01)
	    //    weight[i] = 0.0;
	    //weight[i] = pow(weight[i], 0.9);
	    //weight[i] = 1;
            //double betax=4.0*alphax;
            //double betay=4.0*alphay;
            //double betaxy=16.0*alphaxy;
            //weight[i]=weight[i]-
            //     (exp(-(pow((x[i]-mux),2.0)/betax/betax-
            //                  2.0*betaxy/betax/betax/betay/betay*(x[i]-mux)*(y[i]-muy)+
            //                  pow((y[i]-muy),2.0)/betay/betay)/2.0/(1-pow((betaxy/betax/betay),2.0))))/
            //     (2*PI*betax*betay*sqrt(1-pow((betaxy/betax/betay),2.0)));
            //background=0.0;
	    
            t1 += (x[i]*y[i]*(w[i]-background)*weight[i]);
            t2 += ((w[i]-background)*weight[i]);
            t3 += (x[i]*(w[i]-background)*weight[i]);
            t4 += (y[i]*(w[i]-background)*weight[i]);
            t5 += (x[i]*x[i]*(w[i]-background)*weight[i]);
            t6 += (y[i]*y[i]*(w[i]-background)*weight[i]);
            t7 += (weight[i]*weight[i]);	    
        }
        
	//printf("%lf %lf %lf %lf %lf %lf %lf \n", t1,t2,t3,t4,t5,t6,t7);
        sigmaxy=(t1/t2-t3*t4/t2/t2);
        sigmaxx=(t5/t2-t3*t3/t2/t2);
        sigmayy=(t6/t2-t4*t4/t2/t2);
        mux=(t3/t2);
        muy=(t4/t2);
        signal=t2/t7;
        if (med != 0)
            background=(total-signal)/sizex/sizey;
       //  total = signal + sizex*sizey*background
        rms=sigmaxx+sigmayy;
        if (rms < 0.15) rms=0.15;
        rms=sqrt(rms);
        e1=(sigmaxx-sigmayy)/(sigmaxx+sigmayy);
        e2=(2.0*sigmaxy)/(sigmaxx+sigmayy);
        ellip=sqrt(e1*e1+e2*e2);
        pa=0.5*atan2(e2,e1);
	
        //printf("%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf\n",iteration,signal,mux,muy,rms,e1,e2,background, sigmaxx, sigmayy,sigmaxy, alphax, alphay, alphaxy, pa, ellip);
	//Check if stuck at same value for a while 
        if((iteration > 6) && (abs(sigmaxx - prevsigmaxx) < 1e-3) && (abs(sigmayy - prevsigmayy) < 1e-3) && (abs(sigmaxy - prevsigmaxy) < 1e-3))
        {
	     converged = 1;
	     //printf("ABC \n");
	     break;
        }
	if(iteration %5 == 0)
	{
	    prevsigmaxx = sigmaxx;
	    prevsigmayy = sigmayy;
	    prevsigmaxy = sigmaxy;
	}
	//If actually converged
        if ((abs(alphax-sqrt(2.0*sigmaxx)) < 1e-3) && (abs(alphay-sqrt(2.0*sigmayy)) < 1e-3) && (abs(alphaxy-2.0*sigmaxy) < 1e-3)) {
            converged=1;
            break;
        }
	
	alphax=sigmaxx*2.0;
        if (alphax < 0.15) alphax=0.15;
        alphax=sqrt(alphax);
        alphay=sigmayy*2.0;
        if (alphay < 0.15) alphay=0.15;
        alphay=sqrt(alphay);
        alphaxy=sigmaxy*2.0;
	

    }
    if (converged == 0) {
        //printf("Error: Did not converge.\n");
    } else {
	//This is to fix flux for super small objects(a few pixels)
	if(alphax<0.75 or alphay<0.75)
	{
	    t2=0;
	    t7=0;
	    for (long i = 0; i < count; i++) {
                weight[i]=(exp(-(pow((x[i]-mux),2.0)/alphax/alphax-
                             2.0*alphaxy/alphax/alphax/alphay/alphay*(x[i]-mux)*(y[i]-muy)+
                             pow((y[i]-muy),2.0)/alphay/alphay)/2.0/(1-pow((alphaxy/alphax/alphay),2.0))))/
                (2*PI*alphax*alphay*sqrt(1-pow((alphaxy/alphax/alphay),2.0)));
		t2 += ((w[i]-background)*weight[i]);
		t7 += (weight[i]*weight[i]);
	    }
	    signal=t2/t7;
	}
		
        double backerr = sqrt(signal+background*rms*rms*4*PI)/sqrt(signal);
        printf("\n");
        printf("PhotometricFlux/Counts               %12lf +/- %12lf\n",signal,sqrt(signal)*backerr);
        //printf("AstrometricPositionX/Pixels          %12lf +/- %12lf\n",mux,sqrt(2.0)*rms/sqrt(signal)*backerr);
        //printf("AstrometricPositionY/Pixels          %12lf +/- %12lf\n",muy,sqrt(2.0)*rms/sqrt(signal)*backerr);
        //printf("PsfSize/Pixels                       %12lf +/- %12lf\n",rms,rms/sqrt(signal)*backerr);
        //printf("PsfEllipticityComponent1             %12lf +/- %12lf\n",e1,sqrt(2.0)/sqrt(signal)*backerr);
        //printf("PsfEllipticityComponent2             %12lf +/- %12lf\n",e2,sqrt(2.0)/sqrt(signal)*backerr);
        printf("BackgroundCountRate/(Counts/Pixels)  %12lf +/- %12lf\n",background,sqrt(background*count)/count);
        //printf("aa\n");
	
   	fprintf(fp, "PhotometricFlux/Counts               %12lf  +/- %12lf\n", signal,sqrt(signal)*backerr);
   	fprintf(fp, "AstrometricPositionX/Pixels          %12lf  +/- %12lf\n", mux,sqrt(2.0)*rms/sqrt(signal)*backerr);
	fprintf(fp, "AstrometricPositionY/Pixels          %12lf  +/- %12lf\n",muy,sqrt(2.0)*rms/sqrt(signal)*backerr);
	fprintf(fp, "PsfSize/Pixels                       %12lf  +/- %12lf\n",rms,rms/sqrt(signal)*backerr);
	fprintf(fp, "PsfEllipticityComponent1             %12lf  +/- %12lf\n",e1,sqrt(2.0)/sqrt(signal)*backerr);
	fprintf(fp, "PsfEllipticityComponent2    ,         %12lf  +/- %12lf\n",e2,sqrt(2.0)/sqrt(signal)*backerr);
	fprintf(fp, "BackgroundCountRate/(Counts/Pixels)  %12lf  +/- %12lf\n",background,sqrt(background*count)/count);
   	fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", alphax, alphay, alphaxy, signal, pa, e1, e2, sigmaxx, sigmayy,sigmaxy, mux, muy, rms);
	//fprintf(fp, "%s \n", argv[1]);
    }
    fclose(fp);
}
