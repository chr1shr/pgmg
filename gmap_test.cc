#include <cstdio>

void grid_map(int &lo,int &hi,int ma,bool prd) {
	if(prd) {
		if(lo<-ma||lo>ma-1) lo=-999;
		else lo=lo<0?(lo-(ma&1?2:1))/2:lo>>1;
		if(hi<1||hi>2*ma) hi=-999;
		else hi=(hi+((ma&1)&&hi>ma?3:2))>>1;
	} else {
		if(lo<0||lo>ma-1) lo=-999;
		else lo=(!(ma&1)&&lo==ma-1?lo+1:lo)>>1;
		if(hi<1||hi>ma) hi=-999;
		else hi=(hi+2)>>1;
	}
}

int igrid_map(int &lo,int &hi,int i,int ma,bool prd) {
	if(prd) {
		int j=i*2;
		if(ma&1) {
			int llo,hhi;
			if(j<0) {j+=1;llo=-ma;hhi=0;}
			else if(j>=ma) {j-=1;llo=ma;hhi=2*ma;}
			else {llo=0;hhi=ma;}
			lo=j==llo?llo:j-1;
			hi=j==hhi-1?hhi:j+2;
		} else {
			lo=j-1;
			hi=j+2;
		}
		return j;
	} else {
		if(i<0) {lo=hi=-999;return -999;}
		if(i>=(ma+2)>>1) {lo=hi=-999;return -999;}
		int j=i<<1;
		if(ma&1) {
			lo=j==0?0:j-1;
			hi=j==ma-1?ma:j+2;
		} else {
			if(j==ma) {
				lo=ma-1;hi=ma;return ma-1;
			}
			lo=j==0?0:j-1;
			hi=j==ma-2?ma-1:j+2;
		}
		return j;
	}
}

int main() {
	int c,i,lo,hi;
	for(i=-2;i<=10;i++) {
		lo=i-1;
		hi=i+2;
		grid_map(lo,hi,6,true);
		printf("%3d %3d %3d %3d\n",i-1,i+2,lo,hi);
	}
}
