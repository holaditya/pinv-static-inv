
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct sing_arr {double mat[ 3 ];};
struct matx {double mat[ 3 ][ 3 ];};

struct mat_eqn {double mat[ 4 ];};
struct eig {double mat[ 3 ]};

struct matx trnsps(double mat[ 3 ][ 3 ])
{
    int i, j;
    struct matx mat_trs;
    double trns[ 3 ][ 3 ];

    for(i = 0; i < 3; i++)
    {
   for(j = 0; j < 3; j++)
    {
       mat_trs.mat[ i ][ j ] = mat[ j ][ i ];
    }
    }
    return mat_trs;
}

void print_mat(double mat[3][3])
{
  int i,j;
  for(i = 0; i < 3; i++)
    {
      for(j = 0; j < 3; j++)
      {
        printf( "%lf", mat[i][j]);
        printf(" ");

      }
    printf("\n");
    }
}

struct matx mat_mul(double first[3][3],double second[3][3])
{
  int c,d,k;
  double sum = 0.0;
  struct matx multiply;
  for (c = 0; c < 3; c++) {
    for (d = 0; d < 3; d++) {
      for (k = 0; k < 3; k++) {
        sum = sum + first[c][k]*second[k][d];
  }

  multiply.mat[c][d] = sum;
  sum = 0;
  }
  }

  return multiply;
}

double trace(double mat[3][3])
{
  double sum = 0.0;
  int i;

  for(i = 0; i < 3; i++)
  {
    sum += mat[i][i];
  }
  return sum;
}



struct mat_eqn det_eqn(double a[3][3])
{
  struct matx a2_matx = mat_mul(a, a);
  struct matx a3_matx = mat_mul(a2_matx.mat, a);

  double tr_a = trace(a);
  double tr_a2 = trace(a2_matx.mat);
  double tr_a3 = trace(a3_matx.mat);

  struct mat_eqn eqn;
  eqn.mat[0] = 1;
  eqn.mat[1] = -tr_a;
    //wolfram gets this wrong - does not have the 1/2 factor
    //see: http://math.stackexchange.com/questions/807166/eigenvalues-in-terms-of-trace-and-determinant-for-matrices-larger-than-2-x-2
    //for more information
  eqn.mat[2] = 0.5 * (pow(tr_a, 2) - tr_a2);
    //ditto - wolfram misses the 1/6 factor
  eqn.mat[3] = 1.0 / 6.0 * (pow(tr_a, 3) + (2 * tr_a3) - (3 * tr_a * tr_a2));

  return eqn;
}

struct sing_arr eqn_3(double mat[3][3], double lmd)
{
  printf("\n \n eigen value is %lf\n", lmd);
  int i, j, k, n = 3;
  printf("%lf ", lmd);

  double a[3][4];

    //  Enter the right hand side constants

  for(i = 0; i < n; i++)
  {
    a[i][n] = 0;
  }

    // Enter the coefficients row wise\n
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++)
    {
      a[i][j] = mat[i][j];
      printf("matrix Value (%d) (%d)=%lf\n", i, j, a[i][j]);
    }
  }

  for(i = 0; i < 3; i++)
  {
    a[i][i] = a[i][i] - lmd;
  }

  printf( "matrix Value (%d) (%d)=%lf\n", i, j, a[i][j]);

  struct sing_arr x;
  x.mat[0] = ( a[2][2] * a[1][1] ) - ( a[1][2] * a[2][1] );
  x.mat[1] = ( a[2][2] * a[1][0] )-( a[1][2] * a[2][0] );
  x.mat[1] = ( -1.0 ) * ( x.mat[1] );
  x.mat[2] = ( a[2][1] * a[1][0] )-( a[2][0] * a[1][1] );


  for(i = 0; i < 3; i++)
  printf( " veector=%lf=\n", x.mat[i] );
  printf( "\n" );


  double  sum = sqrt( pow( x.mat[0], 2) + pow( x.mat[1], 2) + pow(x.mat[2], 2));
  printf("\n \n sum %lf \n", sum);
  for(i = 0; i < n; i++)
  {
    x.mat[i] = x.mat[i] / sum;
    printf( "x[%d]= %lf \n", i, x.mat[i] );
  }
 return x;
}



double innner_product(double u[3], double v[3])
{
  double sum = 0.0;
  int i;
  for( i = 0; i < 3; i++ )
    {
      sum=sum + u[i] * v[i];
    }
 return sum;
}

double vector_sign(double u[3], double v[3])
{
 double sum=innner_product( u, v);

 if ( sum > 0 )
{
    return 1.0;
}
  else
  {

    return -1.0;
  }
}

struct sing_arr vector_optimize(double data[ 3 ][ 3 ], double vector[3])//using duke's paper
//https://www.cs.duke.edu/courses/cps111/spring07/notes/12.pdf
// doesn't work
{
  double sum = 0;
  struct sing_arr vec_optimized;
  int chk=1, i;
  double vector_b[3];
  for( i=0; i < 3 ; i++)
  {
    vector_b[ i ]= vector[i];

  }
  for(i = 0; i < 3; i++)
  {
    if(vector_b[i] < 0)
    {
      chk = -1;
      break;
    }
    if(vector_b[i] > 0)
    {
      chk = 1;
      break;
    }
  }

  if(chk == -1)
  {
    for(i = 0; i < 3; i++)
    {
       vector_b[i] = ( -1.0 ) * vector_b[i];
    }
  }

  for( i=0; i < 3 ; i++ )
  {
    vec_optimized. mat[ i ]= vector_b[ i ];
  }
  return vec_optimized;
}




double evaluate_poly(double coeff[4], double point)
{
  double ans = 0;
  int i;
  for( i = 0; i < 4; i++)
  {
    double term = coeff[i] * pow(point, 3 - i);
    ans += term;
  }
return ans;
}

double eval_slope(double coeff[4], double point)
{
  double ans = 0;
  int i;

  for(i = 0; i < 3; i++)
  {
    double term = (3 - i) * coeff[i] * pow(point, 2 - i);
    ans += term;
  }
  return ans;
}

double hack_abs(double x) {
    if (x < 0) { return -x; };
    return x;
}


static const double EPSILON = 0.00000001;
double newton_root(double coeff[4], double begin_x) {
    double current_x = begin_x;
    while(1) {
        double evald_value = evaluate_poly(coeff, current_x);
        double evald_slope = eval_slope(coeff, current_x);
        printf("\nx0: %lf | f(x0): %lf | f'(x0): %lf" , current_x, evald_value, evald_slope);

        if (hack_abs(evald_value) <= EPSILON) {
            return current_x;
        }
        current_x = current_x - evald_value / evald_slope;
    }
};

/*
 double newton_root(double cof[4])
    { int i,j;


        double df[3];

        for(i=0;i<3;i++) // ca
        {   printf(" cof value %lf===",cof[i]);
            df[i]=(3-i)*cof[i];
            printf("dfff %lf \n ",df[i]);
        }
printf("\n");
        double x=300.000;
        double val_dff=0.0;
        double val_fun=0.0;
        for(i=0;i<10;i++)
         { val_fun=0.0;
           val_dff=0.0;
             for(j=0;j<4;j++)
              {

                val_fun=val_fun+cof[j]*pow(x,3-j);

                printf("fun value %lf \n", val_fun );

              }

            printf("fun value %lf \n\n\n\n\n\n", val_fun );
              if (val_fun==0.000)
              {
                  break;
              }
             for(j=0;j<3;j++)
             {  val_dff+=df[j]*pow(x,2-j);

                 printf("diff coff %lf \n diff value %lf\n\n",df[j],val_dff);
             }

        x=x-(val_fun/val_dff);
        printf("\n x is %lf\n",x);

        }

        return x;
    }
*/

void sort_roots( double eigen_roots[3] )
{
 int i = 0, j = 0;
 double temp = eigen_roots[0];
 double tmp_a = 0.0;
 int pos = 0;

 for(i = 0; i < 3; i++)
  {
    temp = eigen_roots[i];
    for( j = i; j < 3; j++)
      {
        if(eigen_roots[ j ] > temp)
           {
                   temp = eigen_roots[ j ];
                   pos = j;
                 }
            }
      tmp_a = eigen_roots[i];
      eigen_roots[ i ] = temp;
      eigen_roots[ pos ] = tmp_a;
    }

 }

struct eig root(double cof[4])
  {
    struct eig eigen;
    eigen.mat [0] = newton_root(cof, 300);

    double qud[3];
    double r2,r3;
    int i;
    double prod;
    double sum = 0;

    if(eigen.mat[0] != 0)
     {
        prod = ( -1 ) * ( cof[3] / cof[0] );
        prod = prod / eigen.mat[0];
        sum=(-1.0 * cof[1] ) / cof[0];
        sum = sum - eigen.mat [0];
        eigen.mat[1] = sum+sqrt( pow( sum, 2 ) - 4 * ( prod ) );
        eigen.mat[1] = eigen.mat[1] / 2;
        eigen.mat[2] = sum-sqrt( pow ( sum, 2 )- 4 * (prod) );
        eigen.mat[2] = eigen.mat[2] / 2;
     }
     else
     {
        prod = cof[2] / cof[0];
        sum = ( -1.0 * cof[1] ) / cof[0];
        sum = sum - eigen.mat[0];
        eigen.mat[1] = sum+sqrt( pow( sum, 2 ) - 4 * ( prod ) );
        eigen.mat[1] = eigen.mat[1] / 2;
        eigen.mat[2] = sum- sqrt( pow( sum, 2)-4 * ( prod ));
        eigen.mat[2] = eigen.mat[2] / 2;

     }


   for(i = 0; i < 3; i++)
     {
       if( hack_abs( eigen.mat[i] ) <= EPSILON)
          eigen.mat[i] = 0.0;
     }
   printf("\n root 1 %lf \n root 2 %lf \n %lf \n n\ ",eigen.mat[1],eigen.mat[2],eigen.mat[0]);
   return eigen;

}


struct sing_arr cross_product( double vec_1[3], double vec_2[3] )
{
    int i;
    struct sing_arr vector;

    vector.mat[0] = vec_1[1] * vec_2[2] - vec_1[2]*vec_2[1];
    vector.mat[1] = ((-1.0) * (vec_1[0] * vec_2[2] - vec_1[2] *vec_2[0]));
    vector.mat[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    printf( "Vec1  \t Vec2  \t Vec3 \n" );
    for(i = 0; i < 3; i++)
     {
       printf( "%lf  \t %lf   \t %lf \n ", vec_1[i], vec_2[i], vector.mat[i] );
     }

    return vector;
};

struct matx svd_u( double mat[3][3], double data[3][3] )
 {
  struct matx mat_trnp = trnsps( mat );
  struct matx mul;
  mul=mat_mul( mat, mat_trnp.mat );
  struct mat_eqn eqn;
  eqn=det_eqn( mul.mat );

  int i,j;
  struct matx u;
  struct sing_arr tmp;

  for(i = 0; i < 4; i++)
   {
     printf( "\n coff in svd_u %lf", eqn.mat[i] );
   }

   struct eig eigen;
   eigen=root(eqn.mat);
   tmp=eqn_3( mul.mat,eigen.mat[0] );
   tmp=vector_optimize(data, tmp.mat );
   double vec_1[3];

   for(j = 0; j < 3; j++)
   {
       vec_1[j] = tmp.mat[j];
   }
   tmp=eqn_3( mul.mat,eigen.mat[1] );
   tmp=vector_optimize(data, tmp.mat);
   double *vec_2 = tmp.mat;
   printf("Vector 1\n");

   for(j = 0; j < 3; j++)
     {
       printf("%lf \t %lf \n ",vec_1[j],tmp.mat[j]);
     }

   for(i = 0; i < 2; i++)
   {
    tmp = eqn_3( mul.mat,eigen.mat[i] );
    tmp=vector_optimize(data, tmp.mat );
    for(j = 0; j < 3; j++)
    {
     u.mat[j][i] = tmp.mat[j];
    }
  }

tmp=cross_product(vec_1,vec_2);
tmp=vector_optimize(data, tmp.mat );

for(j = 0; j < 3; j++)
{
  printf("vec3[%d]=%lf \n ",j,tmp.mat[j]);
  u.mat[j][2] = tmp.mat[j];
  printf("u[%d]=%lf \n ",j,u.mat[j][3]);
}

return u;

}

struct matx svd_sigma(double mat[3][3])
{
  struct matx mat_trnp=trnsps(mat);
  struct matx mul;
  mul = mat_mul(mat, mat_trnp.mat);
  struct mat_eqn eqn;
  eqn = det_eqn(mul.mat);

  int i, j;
  struct matx sigma;
  struct sing_arr tmp;
  struct eig eigen;
  eigen=root(eqn.mat);

  for(i = 0; i < 3; i++ )
  {
    for(j = 0; j < 3; j++ )
    {
      sigma.mat[i][j]=0.0;
    }
  }



  for(i = 0; i < 3; i++ )
  {
    sigma.mat[i][i] = sqrt(eigen.mat[i]);
  }
  return sigma;
}


struct matx eval_pseudo_sigma(double mat[3][3])
{
  struct matx pseudo;
  int i, j;

  printf("\nMatrix insde the pseudo sigma\n");
  print_mat(mat);

  for(i = 0; i < 3; i++)
  {
    for(j = 0;j < 3;j++)
    {
      pseudo.mat[i][j] = 0.0;
    }
  }

  for( i = 0; i < 3 ; i++ )
    {
       if( mat[i][i] != 0 )
        {
           printf( " \n element being processed= %lf \n ", mat[i][i] );
           pseudo.mat[i][i] = 1.0 / mat[i][i];
        }
    }
     return pseudo;

};
struct matx pinv( double mat[3][3], double data[3][3] )
  {
      struct matx u = svd_u( mat, mat );
      struct matx mat_transpose = trnsps( mat );
      struct matx v = svd_u( mat_transpose.mat, data );
      struct matx sigma = svd_sigma( mat );

      printf("U matrix\n");
      print_mat( u.mat );     printf( "\n\n\n\n\n\n" );
      print_mat( v.mat );     printf( "\n\n\n\n\n\n" );
      print_mat( sigma.mat ); printf( "\n\n\n\n\n\n" );


      struct matx pseudo_sigma = eval_pseudo_sigma( sigma.mat );
      printf(" Pseudo S Matrix \n ");
      print_mat( pseudo_sigma.mat );
      struct matx u_transpose = trnsps( u.mat );
      printf("\n U transpose \n");
      print_mat( u_transpose.mat );

      struct matx temp_mul = mat_mul( v.mat, pseudo_sigma.mat );
      printf( "\n Intermediate value \n" );
      print_mat( temp_mul.mat );
      struct matx pseudo_mat = mat_mul( temp_mul.mat, u_transpose.mat );


      printf( " \n\n Value of temp \n " );
      print_mat( temp_mul.mat );
      printf( " \n\n Value of temp \n " );
      print_mat( u_transpose.mat );


      printf( "\n\n\n" );
      print_mat( pseudo_mat.mat );

      return pseudo_mat;
}

int main()
{
  double mat[3][3];
  int i,j;

  for(i= 0; i < 3; i++)
  {
    for(j = 0; j < 3; j++)
    {
      scanf( " %lf ", &mat[ i ][ j ] );
    }
  }

  struct matx pseudo_mat = pinv( mat, mat );
  struct matx mat_check = pinv( pseudo_mat.mat, pseudo_mat.mat);
  return 0;
}
