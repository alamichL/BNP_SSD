library(Rcpp)

cppFunction('NumericVector MvInv_cpp(double eps, double u = 0.5, double alpha = 1, double beta = 1, double gama = 1/2, int N = 3001) {

            NumericVector rg(N);
            rg =   Range(0,N-1) ;
            NumericVector x =   -log( rg * (exp(-10) - exp(-1e-05))/(N-1.0) + exp(-1e-05) );
            
            NumericVector  f =  alpha/Rf_gammafn(1. - gama) * pow(x , -(1. + gama)) * exp(-(u + beta) * x);
            NumericVector  dx = diff(x);
            NumericVector  h = (tail(f,N-1) + head(f,N-1))/2;
            NumericVector  Mv(N);
            
            for (int i = N-2; i >= 0; i--) Mv(i) = Mv(i + 1) + dx(i) * h(i);
            
            double err = 1;
            double w = 0;
            NumericVector v;
            
            while (err > eps){
            w = w + R::rgamma(1,1);
            
            LogicalVector bb = Mv > w;
            int nd=std::find (bb.begin(), bb.end(), 0)-bb.begin(); 
            
            
            v.push_back( x( nd ) ); 
            
            err = min(v)/sum(v); 
            }
            
            return v;
            }')