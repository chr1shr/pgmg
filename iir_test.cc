#include <cstdio>
#include <cstdlib>

void index_range(int *os,int *ws,int vp,int a,int b,int &l,int &u) {
	*ws=*os;
	l=0;
	while(a>=ws[l]) {
		ws[l+1]=ws[l]+os[l+1];
		l++;
		if(l==vp) {
			fputs("Error in transfer indexing [l]\n",stderr);
			exit(1);
		}
	}
	u=l;
	while(b>ws[u]) {
		if(u+1==vp) {
			fputs("Error in transfer indexing [u]\n",stderr);
			exit(1);
		}
		ws[u+1]=ws[u]+os[u+1];
		u++;
	}
	u++;
}

int main() {
	const int vp=4;
	int os[vp]={7,5,4,7},ws[vp],l,u;

	index_range(os,ws,vp,0,23,l,u);

	printf("l=%d u=%d\n",l,u);
	for(int i=0;i<vp;i++) printf("ws[%d]=%d\n",i,ws[i]);
}
