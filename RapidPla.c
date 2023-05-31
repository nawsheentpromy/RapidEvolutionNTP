#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592654
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS (1.2E-07)
#define RNMX (1.0-EPS)

double bnldev(double pp, long n, long *idum); double gammln(double xx); double ran1(long *idum); double randn (double mu, double sigma);

int main(void){

long N1, N2, TIME, REP, interval; double NMU, MIGS, sel, REC, pla;

double epssm, sigma; //random selection values

FILE * par,  * out1; par = fopen("pars1", "r"); out1 = fopen("out2", "w");

fscanf(par, "%ld", &N1);
fscanf(par, "%ld", &N2);
fscanf(par, "%ld", &TIME);
fscanf(par, "%ld", &REP);
fscanf(par, "%ld", &interval);
fscanf(par, "%lf", &NMU);
fscanf(par, "%lf", &MIGS);
fscanf(par, "%lf", &sel);
fscanf(par, "%lf", &REC);
fscanf(par, "%lf", &pla);
fscanf(par, "%lf", &sigma);

double pla1a = 0.0, pla1d = 1.0, pla2a = 1.0, pla2d = 0.0;
double sel1 = sel, sel2 = -sel;

long seed = -1234;

double D1, D2;  // deterministic variables, recombination level //
int one, two;  //  target/modified locus selection variable, mutation level //
long R1, R2;  // reproduction variables for population 1 //
double r1, r2; // reproduction variables for population 1 //
long S1, S2;  // reproduction variables for population 2 //
double s1, s2; // reproduction variables for population 2 //
long k; // time loop //

// time related variables //

long start_time_target = 0, start_time_plastic = 0;
long fixation_time_target = 0, fixation_time_plastic = 0, number_of_fixed_target = 0, number_of_fixed_plastic = 0;
long loss_time_target = 0, loss_time_plastic = 0, number_of_loss_target = 0, number_of_loss_plastic = 0;
double avg_fix_time_target = 0.0, avg_fix_time_plastic = 0.0, avg_loss_time_target = 0.0, avg_loss_time_plastic = 0.0;
double prob_fix_time_target = 0.0, prob_fix_time_plastic = 0.0;

//SELECTION variables //

double avefit1, avefit2, select = 0.0 ;
double a1_m, a1_M, d1_m, d1_M;  // variables for population 1 //
double a2_m, a2_M, d2_m, d2_M;  // variables for popultion 2 //
double fit_a1_m, fit_a1_M, fit_d1_m; // fitness for population 1 //
double fit_a2_m, fit_a2_M, fit_d2_m; // fitness for population 2 //

// MIGRATION variables //
double newN1, newN2, ranpick1, ranpick2, m, m1, m2;
double a1_m_migs=0.0, a1_m_stay, a1_M_migs=0.0, a1_M_stay; //pop1//
double d1_m_migs=0.0, d1_m_stay, d1_M_migs=0.0, d1_M_stay;
double a2_m_migs=0.0, a2_m_stay, a2_M_migs=0.0, a2_M_stay; //pop2//
double d2_m_migs=0.0, d2_m_stay, d2_M_migs=0.0, d2_M_stay;
double a1_mnew, d1_mnew, a1_Mnew, d1_Mnew, a2_mnew, d2_mnew, a2_Mnew, d2_Mnew;

// Cumulative heterozygosity (HL) variables //
double hl_target_locus = 0.0, hl_plastic_locus = 0.0;

long cc, i;

double a,b, ax, bx; //summary variable //

double pop1 = (double)(N1)/(double)(N1+N2);
double pop2 = (double)(N2)/(double)(N1+N2);

int xx; //migration loop

long d; // interval loop

for (i = 0 ; i < REP; i++)
{
  a1_m = 1.000; a2_m = 1.000;
  a1_M = 0.00; d1_m = 0.00; d1_M = 0.00;
  a2_M = 0.00; d2_m = 0.00; d2_M = 0.00;

   one = 0; two = 0;

	  if(ran1(&seed) <= 0.5)
	  {
	  	if(ran1(&seed) < pop1)
	      {
	    	  d1_m = 1.00/(double)N1;
	    	  a1_m = 1.00 - d1_m;
	    	}
	  	else
	      {
	    	  d2_m = 1.00/(double)N2;
	    	  a2_m = 1.00 - d2_m;
	  	  }
	    	one = 1;
	  }

	  else
	  {
	    if(ran1(&seed) < pop1)
	      {
	      	a1_M = 1.00/(double)N1;
	      	a1_m = 1.00-a1_M;
	  	  }
	  	else
	     {
	      	a2_M = 1.00/(double)N2;
	      	a2_m = 1.00-a2_M;
	  	 }
	    	two = 1;
	  }
	  start_time_target = 0;
	  start_time_plastic = 0;

  for( k = 0 ; k< TIME ; k++)
  {
  	////////////////  MUTATION ////////////////////

/////////////  target locus  /////////////

   if (ran1(&seed) < NMU && one == 0)
     {


       if(ran1(&seed) < pop1)
         {
             if(ran1(&seed) < a1_m)
             {
               a1_m = a1_m - 1.00/(double)N1;
               d1_m = d1_m + 1.00/(double)N1;
             }
             else
             {
               a1_M = a1_M - 1.00/(double)N1;
               d1_M = d1_M + 1.00/(double)N1;
             }

         }
       else
         {
             if(ran1(&seed) < a2_m)
             {
               a2_m = a2_m - 1.00/(double)N2;
               d2_m = d2_m + 1.00/(double)N2;
             }

             else
             {
               a2_M = a2_M - 1.00/(double)N2;
               d2_M = d2_M + 1.00/(double)N2;
             }

         }
       one = 1;
       start_time_target = k;
       }


///////////////  plastic locus  /////////////

			   if (ran1(&seed) < NMU && two == 0)
			      {


			        if (ran1(&seed) < pop1)
			           {
			             if(ran1(&seed) < a1_m)
			             {
			               a1_m = a1_m - 1.00/(double)N1;
			               a1_M = a1_M + 1.00/(double)N1;
			             }
			             else
			             {
			               d1_m = d1_m - 1.00/(double)N1;
			               d1_M = d1_M + 1.00/(double)N1;

			             }

			           }
			        else
			           {
			             if(ran1(&seed) < a2_m)
			             {
			               a2_m = a2_m - 1.00/(double)N2;
			               a2_M = a2_M + 1.00/(double)N2;
			             }
			             else
			             {
			               d2_m = d2_m - 1.00/(double)N2;
			               d2_M = d2_M + 1.00/(double)N2;
			             }

			          }
          two = 1;
			    start_time_plastic = k;
        }

///////////////// END MUTATION ////////////////////

  	/////////////////  MIGRATION ////////////////////

  	// if MIGRATION is less than 100, then we will pick individual frequency and check //

    // periodic MIGRATION will take place depending on interval value//

    d = k % interval;

    if (d == 0)
    {

  	if (MIGS < 100.0)
    {
        newN1 = (double)N1;
        newN2 = (double)N2;

        a1_m_stay = a1_m * newN1 ;
        a1_M_stay = a1_M * newN1 ;
        d1_m_stay = d1_m * newN1 ;
        d1_M_stay = d1_M * newN1 ;

        a2_m_stay = a2_m * newN2 ;
        a2_M_stay = a2_M * newN2 ;
        d2_m_stay = d2_m * newN2 ;
        d2_M_stay = d2_M * newN2 ;

        a1_m_migs = 0.0; a1_M_migs = 0.0; d1_m_migs = 0.0; d1_M_migs = 0.0;
        a2_m_migs = 0.0; a2_M_migs = 0.0; d2_m_migs = 0.0; d2_M_migs = 0.0;

      for (xx = 0; xx< MIGS; xx++)
      {
          ranpick1 = ran1(&seed);
          ranpick2 = ran1(&seed);

          a1_m = a1_m_stay/newN1 ;
          a1_M = a1_M_stay/newN1 ;
          d1_m = d1_m_stay/newN1 ;

          a2_m = a2_m_stay/newN2 ;
          a2_M = a2_M_stay/newN2 ;
          d2_m = d2_m_stay/newN2 ;

          if (ranpick1 <= a1_m)
          {
            a1_m_migs = a1_m_migs + 1.0;
            a1_m_stay = a1_m_stay - 1.0;
          }

          else if (ranpick1 > a1_m && ranpick1 <= (a1_m + a1_M))
          {
            a1_M_migs = a1_M_migs + 1.0;
            a1_M_stay = a1_M_stay - 1.0;
          }

          else if (ranpick1 > (a1_m + a1_M) && ranpick1 <= (a1_m + a1_M + d1_m))
          {
            d1_m_migs = d1_m_migs + 1.0;
            d1_m_stay = d1_m_stay - 1.0;
          }

          else
          {
            d1_M_migs = d1_M_migs + 1.0;
            d1_M_stay = d1_M_stay - 1.0;
          }

          if (ranpick2 <= a2_m)
          {
            a2_m_migs = a2_m_migs + 1.0;
            a2_m_stay = a2_m_stay - 1.0;
          }

          else if (ranpick2 > a2_m && ranpick2 <= (a2_m + a2_M))
          {
            a2_M_migs = a2_M_migs + 1.0;
            a2_M_stay = a2_M_stay - 1.0;
          }

          else if (ranpick2 > (a2_m + a2_M) && ranpick2 <= (a2_m + a2_M + d2_m))
          {
            d2_m_migs = d2_m_migs + 1.0;
            d2_m_stay = d2_m_stay - 1.0;
          }

          else
          {
            d2_M_migs = d2_M_migs + 1.0;
            d2_M_stay = d2_M_stay - 1.0;
          }

 	  newN1 = newN1 - 1.0;
          newN2 = newN2 - 1.0;
      
} // end for loop under migs<100

      a1_mnew = (a1_m_stay + a2_m_migs)/(double)N1;
      a1_Mnew = (a1_M_stay + a2_M_migs)/(double)N1;
      d1_mnew = (d1_m_stay + d2_m_migs)/(double)N1;
      d1_Mnew = (d1_M_stay + d2_M_migs)/(double)N1;

      a2_mnew = (a2_m_stay + a1_m_migs)/(double)N2;
      a2_Mnew = (a2_M_stay + a1_M_migs)/(double)N2;
      d2_mnew = (d2_m_stay + d1_m_migs)/(double)N2;
      d2_Mnew = (d2_M_stay + d1_M_migs)/(double)N2;

      a1_m = a1_mnew;
      d1_m = d1_mnew;
      a1_M = a1_Mnew;
      d1_M = d1_Mnew;
      a2_m = a2_mnew;
      d2_m = d2_mnew;
      a2_M = a2_Mnew;
      d2_M = d2_Mnew;

    } //end if migs<100


    else if (MIGS >= 100.0)
    {

      //m = MIGS/(double)(N1+N2);
      m1 = MIGS/(double)N1;
      m2 = MIGS/(double)N2;

      a1_mnew = a1_m*(1.0-m1)+a2_m*m1;
      d1_mnew = d1_m*(1.0-m1)+d2_m*m1;

      a1_Mnew = a1_M*(1.0-m1)+a2_M*m1;
      d1_Mnew = d1_M*(1.0-m1)+d2_M*m1;

      a2_mnew = a2_m*(1.0-m2)+a1_m*m2;
      d2_mnew = d2_m*(1.0-m2)+d1_m*m2;

      a2_Mnew = a2_M*(1.0-m2)+a1_M*m2;
      d2_Mnew = d2_M*(1.0-m2)+d1_M*m2;

    	a1_m = a1_mnew;
    	d1_m = d1_mnew;
    	a1_M = a1_Mnew;
    	d1_M = d1_Mnew;
    	a2_m = a2_mnew;
    	d2_m = d2_mnew;
    	a2_M = a2_Mnew;
    	d2_M = d2_Mnew;

    }   // end else if

    else
    {
      printf("Migration is not working");
    }

  }


			///////////////////   SELECTION     ///////////////////
        {
  				epssm = sigma*randn(0.0, 1.0);
        }
				
			//////////////////    population 1     //////////////////

			avefit1 = a1_m*(1.0 +(1.0+epssm)*sel1) + d1_m*(1.0 - (1.0+epssm)*sel1) + a1_M*(1.0 +(1.0+epssm)*sel1*(1.0-pla1a)) +  d1_M*(1.0 - (1.0+epssm)*sel1 *(1.0-pla1d));
			
			fit_a1_m= (1.0 + (1.0+epssm)*sel1)/avefit1;
			fit_d1_m= (1.0 - (1.0+epssm)*sel1)/avefit1;
			fit_a1_M= (1.0 + (1.0+epssm)*sel1*(1.0-pla1a))/avefit1;
			
			if(a1_m > 0.0) a1_m = fit_a1_m*a1_m;
			if(a1_M > 0.0) a1_M = fit_a1_M*a1_M;
			if(d1_m > 0.0) d1_m = fit_d1_m*d1_m;
			if(d1_M > 0.0) d1_M = 1.0 - a1_m - a1_M - d1_m;

			////////////////    population 2     //////////////////

			avefit2 = a2_m*(1.0 + (1.0+epssm)*sel2) + d2_m*(1.0 - (1.0+epssm)*sel2) + a2_M*(1.0 + (1.0+epssm)*sel2*(1.0-pla2a)) +  d2_M*(1.0 - (1.0+epssm)*sel2 *(1.0-pla2d));
			
			fit_a2_m= (1.0 + (1.0+epssm)*sel2)/avefit2;
			fit_d2_m= (1.0 - (1.0+epssm)*sel2)/avefit2;
			fit_a2_M= (1.0 + (1.0+epssm)*sel2*(1.0-pla2a))/avefit2;
			
			if(a2_m > 0.0) a2_m = fit_a2_m*a2_m;
			if(a2_M > 0.0) a2_M = fit_a2_M*a2_M;
			if(d2_m > 0.0) d2_m = fit_d2_m*d2_m;
			if(d2_M > 0.0) d2_M = 1.0 - a2_m - a2_M - d2_m;

			///////////////////   END SELECTION     ///////////////////

    /////////////  RECOMBINATION   ///////////////

	/////////////  population 1  ///////////////

	if((a1_m+a1_M)*(d1_m+d1_M)*(a1_m+d1_m)*(a1_M+d1_M) > 0.0)
	{

	  D1 = ((a1_m*d1_M)-(a1_M*d1_m))*REC;

	  a1_m = a1_m - D1;
	  a1_M = a1_M + D1;
	  d1_m = d1_m + D1;
	  d1_M = d1_M - D1;

	}

	/////////////  population 2  ///////////////

	if((a2_m+a2_M)*(d2_m+d2_M)*(a2_m+d2_m)*(a2_M+d2_M) > 0.0)
	{

	  D2 = ((a2_m*d2_M)-(a2_M*d2_m))*REC;

	  a2_m = a2_m - D2;
	  a2_M = a2_M + D2;
	  d2_m = d2_m + D2;
	  d2_M = d2_M - D2;

	}

	///////////// END RECOMBINATION  ///////////////


	/////////////  REPRODUCTION  ///////////////


    //////////// population 1  //////////////

    if(d1_m > 0.0) r1 = d1_m/(d1_m+d1_M);
    else r1 = 0.0;

    if(a1_M > 0.0) r2 = a1_M/(d1_m+a1_M+d1_M);
    else r2 = 0.0;

    a1_m = bnldev(a1_m, N1, &seed);
    R1 = N1 - (long)a1_m;

    a1_M = bnldev(r2, R1, &seed);
    R2 = R1 - (long)a1_M;

    d1_m = bnldev(r1, R2, &seed);

    d1_M = ((double)N1 - a1_m - a1_M - d1_m)/(double)N1;
    a1_m = a1_m/(double)N1;
    a1_M = a1_M/(double)N1;
    d1_m = d1_m/(double)N1;

    //////////// population 2  //////////////

    if(d2_m > 0.0) s1 = d2_m/(d2_m+d2_M);
    else s1 = 0.0;

    if(a2_M > 0.0) s2 = a2_M/(d2_m+a2_M+d2_M);
    else s2 = 0.0;

    a2_m = bnldev(a2_m, N2, &seed );
    S1 = N2 - (long)a2_m;

    a2_M = bnldev(s2, S1, &seed );
    S2 = S1 - (long)a2_M;

    d2_m = bnldev(s1, S2, &seed );

    d2_M = ((double)N2 - a2_m - a2_M - d2_m)/(double)N2;
    a2_m = a2_m/(double)N2;
    a2_M = a2_M/(double)N2;
    d2_m = d2_m/(double)N2;

    ///////////////// END REPRODUCTION /////////////////


    //////////////// SUMMARIES ////////////////////

	hl_target_locus  += 2*( (a1_m+a1_M)*(1.0-(a1_m+a1_M))*pop1+(a2_m+a2_M)*(1.0-(a2_m+a2_M))*(1.0-pop1));
	hl_plastic_locus += 2*( (a1_m+d1_m)*(1.0-(a1_m+d1_m))*pop1+(a2_m+d2_m)*(1.0-(a2_m+d2_m))*(1.0-pop1));

	ax = (((a1_m+a1_M)*(double)N1 + (a2_m+a2_M)*(double)N2))/((double)(N1+N2));
	bx = (((a1_m+d1_m)*(double)N1 + (a2_m+d2_m)*(double)N2))/((double)(N1+N2));


	//  target  //

	if((ax*(1-ax) == 0.0) && one == 1)
	  {

	  if(ax == 0.0)
	    {
	      fixation_time_target += k-start_time_target;
	      number_of_fixed_target +=1;
	    }

	  else if(ax == 1.0)
	    {
	      loss_time_target += k-start_time_target;
	      number_of_loss_target +=1;
	    }
	  else
	    {
	      printf("Target Not Working\n");
	    }
	  one = 2;

	  }


	  //  plastic  //

	    if((bx*(1-bx) == 0.0) && two == 1)
	  {

	  if(bx == 0.0)
	    {
	      fixation_time_plastic += k-start_time_plastic;
	      number_of_fixed_plastic +=1;
	    }

	  else if(bx == 1.0)
	    {
	      loss_time_plastic += k-start_time_plastic;
	      number_of_loss_plastic +=1;
	    }

	  else
	    {
	      printf("Plastic Not Working\n");
	    }
	  two = 2;

	  }

	  if(one > 1 && two > 1)
	  {
	  break;
	  }


  } // end TIME


	if( i==99999 && hl_target_locus/100000.0 > 20.0 && hl_plastic_locus/100000.0 > 20.0)
	{
	 REP = 100000;
	 break;
	}

} // end REP

	avg_fix_time_target = (double)fixation_time_target/(double)number_of_fixed_target;
	avg_fix_time_plastic = (double)fixation_time_plastic/(double)number_of_fixed_plastic;

	avg_loss_time_target = (double)loss_time_target/(double)number_of_loss_target;
	avg_loss_time_plastic = (double)loss_time_plastic/(double)number_of_loss_plastic;

	prob_fix_time_target = (double)number_of_fixed_target/(double)(REP);
	prob_fix_time_plastic = (double)number_of_fixed_plastic/(double)(REP);


fprintf(out1, "N1 %ld N2 %ld MIGS %.0f sel1 %.3f sel2 %.3f pla1a %.3f pla1d %.3f pla2a %.3f pla2d %.3f  hl_target_locus %.4f hl_plastic_locus %.4f avg_fix_time_target %.4f avg_fix_time_plastic %.4f avg_loss_time_target %.4f avg_loss_time_plastic %.4f prob_fix_target  %.8f prob_fix_plastic = %.8f\n", N1, N2, MIGS, sel1, sel2, pla1a, pla1d, pla2a, pla2d, hl_target_locus/(double)REP ,hl_plastic_locus/(double)REP, avg_fix_time_target, avg_fix_time_plastic, avg_loss_time_target, avg_loss_time_plastic, prob_fix_time_target, prob_fix_time_plastic);
fclose(par);
fclose(out1);

return 0;

} // END MAIN //

double randn (double mu, double sigma){
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}


double bnldev(double pp, long n, long *idum)
{

long j;

static long nold =(-1);

double am, em, g, angle, p, bnl, sq, t, y;
static double pold=(-1.0), pc, plog, pclog, en, oldg;

p=(pp <= 0.5 ? pp : 1.0-pp);

am = n*p;

if(n<25){
	bnl = 0.0;

	for(j=1; j<=n; j++)
		if(ran1(idum) < p) ++bnl;
}else if (am < 1.0) {
  	g=exp(-am);
  	t = 1.0;
	for (j=0; j<=n; j++){
		t *= ran1(idum);
		if (t < g) break;
}
bnl=(j<=n ? j:n);
}else {


if (n != nold) {
	en = n;
	oldg = gammln(en+1.0);
	nold = n;
}

if (p != pold) {
	pc = 1.0-p;
	plog = log(p);
	pclog = log(pc);
	pold = p;
}

sq = sqrt(2.0*am*pc);
do{
	do{
		angle = PI * ran1(idum);
		y = tan(angle);
		em = sq*y+am;
	} while (em < 0.0 || em >= (en + 1.0));

	em = floor(em);
	t = 1.2 * sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
}while (ran1(idum) > t);


bnl = em;

}

if (p != pp) bnl = n - bnl;
return bnl;

}

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	long j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double ran1(long *idum)
{
	long j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

