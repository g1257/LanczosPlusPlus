#include "GitRevision.h"

int main(int argc,char *argv[])
{
	PsimagLite::GitRevision gitrev("./","lanczos");
	std::cout<<gitrev;
	PsimagLite::GitRevision gitrev2("../../PsimagLite/","psimagLite");
	std::cout<<gitrev2;

}
