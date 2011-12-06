#if !defined _grv94_h_
#define _grv94_h_

double ss (double Q2)//to jest poczatek GRV94
{
  double a2=232.2*232.2;
  double mi2=230000;
  if (Q2<=232000)//wszystko jest w MeV...
    return 0;
  else
    return log(log(Q2/a2)/log(mi2/a2));//log naturalny
}



double uv (double xx, double Q2)//u(val), dobre unormowanie!
{
  double s=ss(Q2);
  double s2=s*s;

//do wzoru na u(val)
double a = 0.59-0.024*s;
double b = 0.131+0.063*s;

double N = 2.284+0.802*s+0.055*s2;

double A = -0.449-0.138*s-0.076*s2;

double B = 0.213+2.669*s-0.728*s2;

double C = 8.854-9.135*s+1.979*s2;

double D = 2.997+0.753*s-0.076*s2;

return N*pow(xx,a)*(1+A*pow(xx,b)+xx*(B+C*sqrt(xx)))*pow(1-xx,D);

}



double dv (double xx, double Q2)//d(val), dobre unormowanie!
{ 
  double s=ss(Q2);
  double s2=s*s;
  
  //do wzoru na d(val)

  
  double a = 0.376;
  double b = 0.486+0.062*s;

  double N = 0.371+0.083*s+0.039*s2;

  double A = -0.509+3.31*s-1.248*s2;

  double B = 12.41-10.52*s+2.267*s2;

  double C = 6.373-6.208*s+1.418*s2;
  
  double D = 3.691+0.799*s-0.071*s2;

  
  return N*pow(xx,a)*(1+A*pow(xx,b)+ xx*(B+C*sqrt(xx)))*pow(1-xx,D);
}



double del (double xx, double Q2)//dbar-ubar
{
  double s=ss(Q2);
  double s2=s*s;
  
  //do wzoru na dbar-ubar
  
  
  double a = 0.409-0.005*s;
  double b = 0.799+0.071*s;
  
  double N = 0.082+0.014*s+0.008*s2;
  
  double A = -38.07+36.13*s-0.656*s2;
  
  double B = 90.31-74.15*s+7.645*s2;
  
  double C = 0;
  
  double D = 7.486+1.217*s-0.159*s2;
  
  
  return N*pow(xx,a)*(1+A*pow(xx,b)+xx*B)*pow(1-xx,D);
}



double sb (double xx, double Q2)//s=s(bar)
{
  double s=ss(Q2);
  double s2=s*s; 

//do wzoru na s=s(bar)

double alfa = 0.914;
double beta = 0.577;

double a = 1.798-0.596*s;

double A = -5.548+3.669*sqrt(s)-0.616*s;

double B = 18.92-16.73*sqrt(s)+5.168*s;

double D = 6.379-0.35*s+0.142*s2;

double E = 3.981+1.638*s;

double Eprim = 6.402;

return pow(s,alfa)/pow(log(1.0/xx),a)*(1+A*sqrt(xx)+B*xx )*
pow(1-xx,D)*exp(-E+sqrt(Eprim*pow(s,beta)*log(1.0/xx)  ) );
}


double udb (double xx, double Q2)//ubar+dbar
{

  double s=ss(Q2);
  double s2=s*s;

//do ubar+dbar

double alfa = 1.451;
double beta = 0.271;

double a = 0.41-0.232*s;
double b = 0.534-0.457*s;

double A = 0.89-0.14*s;

double B = -0.981;

double C = 0.32+0.683*s;

double D = 4.752+1.164*s+0.286*s2;

double E = 4.119+1.713*s;

double Eprim = 0.682+2.978*s;


return (pow(xx,a)*(A+B*xx+C*xx*xx)*pow(log(1.0/xx),b)+
pow(s,alfa)*exp(-E+sqrt(Eprim*pow(s,beta)*log(1.0/xx))))*pow(1-xx,D);
}


double dbar (double xx, double Q2)
{return (del(xx, Q2)+udb(xx,Q2))/2;}
double ubar (double xx, double Q2)
{return (-del(xx, Q2)+udb(xx,Q2))/2;}


//double F3_n (double xx, double Q2)
//{return 2*(uv(xx,Q2)+sb(xx,Q2)+ubar(xx,Q2)-dbar(xx,Q2));}//wyrazenia na funkcje struktury
//double F3_p (double xx, double Q2)
//{return 2*(dv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)-ubar(xx,Q2));}
//double F2_n (double xx, double Q2)
//{return 2*xx*(uv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)+ubar(xx,Q2));}
//double F2_p (double xx, double Q2)
//{return 2*xx*(dv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)+ubar(xx,Q2));}
//double F1_n (double xx, double Q2)
//{return F2_n(xx,Q2)/2/xx;}
//double F1_p (double xx, double Q2)
//{return F2_p(xx,Q2)/2/xx;}


//zmiana KG

double F3_cc_nu_n  (double xx, double Q2) {return 2*(uv(xx,Q2)+sb(xx,Q2)+ubar(xx,Q2)-dbar(xx,Q2))/xx;}//wyrazenia na funkcje struktury
double F2_cc_nu_n  (double xx, double Q2) {return 2*(uv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)+ubar(xx,Q2));}
double F1_cc_nu_n  (double xx, double Q2) {return F2_cc_nu_n(xx,Q2)/2/xx;}
double F3_cc_anu_n (double xx, double Q2) {return 2*(dv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)-ubar(xx,Q2))/xx;}//wyrazenia na funkcje struktury
double F2_cc_anu_n (double xx, double Q2) {return 2*(dv(xx,Q2)+sb(xx,Q2)+ubar(xx,Q2)+dbar(xx,Q2));}
double F1_cc_anu_n (double xx, double Q2) {return F2_cc_nu_n(xx,Q2)/2/xx;}

double F3_cc_nu_p  (double xx, double Q2) {return 2*(dv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)-ubar(xx,Q2))/xx;}
double F2_cc_nu_p  (double xx, double Q2) {return 2*(dv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)+ubar(xx,Q2));}
double F1_cc_nu_p  (double xx, double Q2) {return F2_cc_nu_p(xx,Q2)/2/xx;}
double F3_cc_anu_p (double xx, double Q2) {return 2*(uv(xx,Q2)+sb(xx,Q2)+ubar(xx,Q2)-dbar(xx,Q2))/xx;}
double F2_cc_anu_p (double xx, double Q2) {return 2*(uv(xx,Q2)+sb(xx,Q2)+ubar(xx,Q2)+dbar(xx,Q2));}
double F1_cc_anu_p (double xx, double Q2) {return F2_cc_nu_p(xx,Q2)/2/xx;}


const double g_V_u = 1/2. - 4/3.*sin_2_theta_W; 
const double g_V_d = 2/3.*sin_2_theta_W -1/2.;
const double g_A_u =  1/2.; 
const double g_A_d = -1/2.;


//funkcje struktury dla NC by JN
double F3_nc_n(double xx, double Q2) 
{
return 2*(uv(xx,Q2)+sb(xx,Q2)+ubar(xx,Q2)-dbar(xx,Q2))/xx;
}
double F2_nc_n(double xx, double Q2) 
{
return 2*(uv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)+ubar(xx,Q2));
}
double F1_nc_n(double xx, double Q2) 
{
return F2_cc_nu_n(xx,Q2)/2/xx;
}

double F3_nc_p(double xx, double Q2) 
{
return 2*(dv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)-ubar(xx,Q2))/xx;
}
double F2_nc_p(double xx, double Q2) 
{
return 2*(dv(xx,Q2)+sb(xx,Q2)+dbar(xx,Q2)+ubar(xx,Q2));
}
double F1_nc_p(double xx, double Q2) 
{
return F2_cc_nu_p(xx,Q2)/2/xx;
}


#endif
