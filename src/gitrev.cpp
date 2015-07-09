#include "GitRevision.h"

int main(int,char **)
{
	PsimagLite::GitRevision gitrev("./","lanczos");
	std::cout<<gitrev;
	PsimagLite::GitRevision gitrev2("../../PsimagLite/","psimagLite");
	std::cout<<gitrev2;

}
