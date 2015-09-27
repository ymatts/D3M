//' @useDynLib D3M
//' @importFrom Rcpp sourceCpp


#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
NumericVector quantileCpp(NumericVector x, NumericVector probs){
  
  int xn = x.size();
  
  int pn = probs.size();
  
  double bin = (double)(1.0 / xn);
  
  NumericVector xs = clone(x);
  
  std::sort(xs.begin(), xs.end());
  
  IntegerVector pos(pn);
  
  int temp = 0;
  for(int j = 0; j < pn; j++){
    
    temp = floor(xn * probs[j]) - 1;
    if(temp < 0){
      temp = 0;
    }
    pos[j] = temp;
    //printf("pos[%d] = %d\n",j,pos[j]);
    
  }
  
  NumericVector out(pn);
  
  //Initializing output vector
  
  
  for(int i = 0; i < pn; i++){
    out[i] = 0;
  }  
  
  int index = 0;
  
  double accum_bin = 0;
  
  
  for(int j = 0; j < pn; j++){
    
    if(probs[j] == 0.0){
      
      out[j] = min(xs);
      
      continue;
    }
    
    if(probs[j] == 1.0){
      
      out[j] = max(xs);
      
      continue;
    }
    
    index = pos[j];
    
    accum_bin = bin * (index + 1);
    
    //printf("xs %f index %d probs %f accum_bin %f\n",xs[index],index, probs[j], accum_bin);
    out[j] = xs[index] + ((probs[j] - accum_bin) / bin) * (xs[index + 1] - xs[index]);
    
  }
  
  return out;
  //return result;
}
// [[Rcpp::export]]
double wasserCpp(NumericVector x, NumericVector y, int paranum = 101, int q = 2){
  //int len = x.size();
  //printf("paranum=%d q=%d\n",paranum,q);
  
  //NumericVector iqr_x = quantileCpp(x,IQR);
  //NumericVector iqr_y = quantileCpp(y,IQR);
  
  //NumericVector x2 = x[(iqr_x[0] <= x) & (x <= iqr_x[1])];
  //NumericVector y2 = y[(iqr_y[0] <= y) & (y <= iqr_y[1])];
  
  NumericVector pseq(paranum);
  NumericVector xq(paranum), yq(paranum);
  
  for(int i = 0; i < paranum; i++){
    pseq[i] = (double)i / (paranum - 1);
    //printf("%.3f\t",pseq[i]);
  }
  //printf("\n");
  
  //x = x[!is_na(x)];
  //y = y[!is_na(y)];
  
  xq = quantileCpp(x, pseq);
  yq = quantileCpp(y, pseq);
  
  double d = 0;
  
  for(int i = 0; i < paranum; i++){
    //printf("%.3f\t",xq[i]);
    //printf("%3.f\t",yq[i]);
    d += pow(xq[i] - yq[i],q);
  }
  return d;
}

// [[Rcpp::export]]
NumericVector wasserCpp_mat(NumericMatrix xMat, NumericMatrix yMat, int paranum = 101, int q = 2){
  
  int nr_xMat = xMat.nrow();
  int nc_xMat = xMat.ncol();
  int nr_yMat = yMat.nrow();
  int nc_yMat = yMat.ncol();
  //printf("x %d %d y %d %d\n",nr_xMat,nc_xMat,nr_yMat,nc_yMat);
  
  if(nr_xMat != nr_yMat){
    
    //printf("Error: size of matrix is different");
    
    double d = NA_REAL;
    
    return d;
  }
  
  int npos = nr_xMat;
  
  //printf("npos %d\n",npos);
  
  NumericVector x(nc_xMat);
  NumericVector y(nc_yMat);
  
  //printf("hoge1\n");
  
  NumericVector d_vec(npos);
    
  for(int pos = 0; pos < npos; pos++){
    
    //printf("hoge2\n");
    
    for(int j = 0; j < nc_xMat; j++){
      //for(int j = 0; j < paranum; j++){
      //printf("hoge3\n");
      x[j] = xMat(pos,j);
      //printf("x%d %.3f\t",j,x[j]);
      
    }
    
    for(int j = 0; j < nc_yMat; j++){
      //for(int j = 0; j < paranum; j++){
      
      y[j] = yMat(pos,j);
      
      //printf("y%d %.3f\t",j,y[j]);
    }
    
    //int len = x.size();
    double d = 0;
    
    NumericVector pseq(paranum);
    NumericVector xq(paranum), yq(paranum);
    
    for(int i = 0; i < paranum; i++){
      
      pseq[i] = (double)i / (paranum - 1);
      
    }
    
    //x = x[!is_na(x)];
    //y = y[!is_na(y)];
    //NumericVector iqr_x = quantileCpp(x,IQR);
    //NumericVector iqr_y = quantileCpp(y,IQR);
    
    //NumericVector x2 = x[(iqr_x[0] <= x) & (x <= iqr_x[1])];
    //NumericVector y2 = y[(iqr_y[0] <= y) & (y <= iqr_y[1])];
    
    
    xq = quantileCpp(x, pseq);
    yq = quantileCpp(y, pseq);
    
    for(int i = 0; i < paranum; i++){
      //printf("%.3f\t",xq[i]);
      d += pow(xq[i] - yq[i],q);
      
    }
    
    d_vec[pos] = d;
    
  }
  
  return d_vec;
}

// [[Rcpp::export]]
List permCpp(NumericMatrix casesMat, NumericMatrix controlMat, NumericMatrix shuffleID, NumericVector d,int bsn = 10000, int qn = 101, int q = 2){
  
  int nr_casesMat = casesMat.nrow();
  
  int nc_casesMat = casesMat.ncol();
  //int nr_controlMat = controlMat.nrow();
  
  int nc_controlMat = controlMat.ncol();
  
  //printf("nr_caseMat %d\nnc_caseMat %d\nnr_controlMat %d\nnc_controlMat %d\n",nr_casesMat,nc_casesMat,nr_controlMat,nc_controlMat);
  
  int npos = nr_casesMat;
  
  NumericVector pval_vec(npos);
  
  NumericVector cases(nc_casesMat);
  
  NumericVector control(nc_controlMat);
  
  //List lbootd(npos);
  
  //List ldata(npos);
  
  //List lcasedata(bsn);
  //List lcontroldata(bsn);
  
  for(int pos = 0; pos < npos; pos++){
    
    NumericVector cases(nc_casesMat);
    NumericVector control(nc_controlMat);
    
    for(int j = 0; j < nc_casesMat; j++){
      //printf("%.3f:",casesMat(pos,j));
      cases[j] = casesMat(pos,j); 
      //printf("%.3f\n",cases[j]);
    }
    
    for(int j = 0; j < nc_controlMat; j++){
      //printf("%.3f:",controlMat(pos,j));
      
      control[j] = controlMat(pos,j);
      
      //printf("%.3f\n",control[j]);
      // printf("i,j %d %d: %f\n",pos,j,control[j]);
    }
    
    //lcases[pos] = cases;
    //lcontrol[pos] = control;
    //cases = cases[!is_na(cases)];
    //control = control[!is_na(control)];
    
    int ncases = cases.size();
    int ncontrol = control.size();
    int nsample = ncases + ncontrol;
    
    //printf("sample size:\ncase %d control %d total %d\n",ncases,ncontrol,nsample);
    
    NumericVector data(nsample);
    
    int cnt = 0;
    
    for(int i = 0; i < ncases; i++){
      
      data[cnt] = cases[i];
      
      cnt++;
    }
    
    for(int i = 0; i < ncontrol; i++){
      data[cnt] = control[i];
      
      cnt++;
    }
    
    // ldata[pos] = data;
    NumericVector index(nsample);
    
    for(int i = 0; i < nsample; i++){
      
      index[i] = i;
    }
    
    NumericVector bootd(bsn);
    
    //NumericVector index2;
    
    
    
    //int nonaN = 0;
    
    //NumericMatrix kk(nsample,bsn);
    
    for(int i = 0; i < bsn; i++){
      
      NumericVector casedata(ncases);
      
      NumericVector controldata(ncontrol);
      
      //NumericVector index2 = randomShuffle(index);
      
      
      int cnt = 0;
      
      for(int j = 0; j < ncases; j++){
        
        int k;
        //k = index2[cnt];
        
        k = shuffleID(cnt,i) - 1;
        
        //kk(cnt,i) = k;
        
        casedata[j] = data[k];
        
        //printf("casedata[%d] = %.3f = data[%d]\n",j,casedata[j],k);
        
        cnt++;
        
      }
      
      //lcasedata[i] = casedata;
      
      for(int j = 0; j < ncontrol; j++){
        
        int k;
        //k = index2[cnt];
        
        k = shuffleID(cnt,i) - 1;
        
        //kk(cnt,i) = k;
        
        controldata[j] = data[k];
        
        //printf("controldata[%d] = %.3f = data[%d]\n",j,controldata[j],k);
        
        cnt++;
        
      }
      
      //lcontroldata[i] = controldata;
      
      double a = wasserCpp(casedata, controldata, qn, q);
      
      //printf("a[%d] = %.3f\n",i,a);
      
      /*
      if(isnan(a) || isinf(a)){  
      
      bootd[i] = NAN;
      
      }else{
      
      bootd[i] = a;
      
      nonaN++;
      }
      */
      bootd[i] = a;
    }
    
    //lshuffleID[pos] = kk;
    
    //printf("nonaN:%d\n",nonaN);
    
    /*
    NumericVector newbootd(nonaN);
    
    cnt = 0;
    
    if(nonaN < bsn){
    for(int l = 0; l < bootd.size(); l++){
    
    if(!isnan(bootd[l]) && !isinf(bootd[l])){
    
    newbootd[cnt] = bootd[l];
    
    cnt++;
    }      
    }
    }else{
    newbootd = clone(bootd);
    }
    */
    //lbootd[pos] = bootd;
    
    //NumericVector temp = newbootd[newbootd >= d[pos]];
    NumericVector temp = bootd[bootd >= d[pos]];
    
    int ln = 0;
    
    //for(int l = 0; l < newbootd.size(); l++){
    for(int l = 0; l < bootd.size(); l++){
      
      //if(newbootd[l] >= d[pos]){
      if(bootd[l] >= d[pos]){  
        
        ln++;
        
      }
    }
    //printf("ln:%d\n",ln);
    
    //double pval = (double)temp.size() / (double)newbootd.size();
    
    double pval = (double) temp.size() / (double) bootd.size();
    
    //printf("pval=temp.size()/newbootd.size()=%.3f/%.3f= %.3f\n",(double)temp.size(), (double)newbootd.size(),pval);
    //printf("pval=temp.size()/newbootd.size()=%.3f/%.3f= %.3f\n",(double)temp.size(), (double)bootd.size(),pval);
    
    
    if(pval < (double) 1 / bsn){
      
      //if(pval == 0.0){  
      //int threshold_ind = (int)bsn * 0.995;
      
      //sort(newbootd.begin(),newbootd.end());
      
      //double threshold = newbootd[threshold_ind];
      
      NumericVector q(1);
      q = 0.995;
      
      //NumericVector thresholdv = quantileCpp(newbootd, 0.995);
      //NumericVector thresholdv = quantileCpp(newbootd, q);
      NumericVector thresholdv = quantileCpp(bootd, q);
      
      double threshold = thresholdv[0];
      
      //printf("threshold %f\n",threshold);
      
      //NumericVector ptemp = newbootd[newbootd > threshold];
      NumericVector ptemp = bootd[bootd > threshold];
      
      //printf("ptemp.size()=%d\n",ptemp.size());
      
      for(int j = 0; j < ptemp.size(); j++){
        
        ptemp[j] -= threshold;
        
      }
      
      
      double lambda = 0;
      
      for(int j =0; j < ptemp.size(); j++){
        
        lambda += ptemp[j];
        
      }
      
      //printf("lambda %f\n",lambda);
      
      lambda = 1 / (lambda / ptemp.size());
      
      //printf("Estimated lambda: %f\n",lambda);
      
      //double param = -1 * lambda * (d[pos] - threshold);
      double param = -1 * lambda * d[pos];
      
      double param2 = -1 * lambda * threshold;
      
      double pval_threshold = exp(param2);
      
      double r = pval_threshold / 0.005;
      
      pval = exp(param) / r;
      
      //printf("semi-parametric p-value (lambda = %.3f) %.3e\n",lambda,pval);
      
    }
    
    pval_vec[pos] = pval;
    
  }
  
  return List::create(pval_vec, d, casesMat, controlMat);
  
}


