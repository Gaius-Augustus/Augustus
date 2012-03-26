// fixed2varialeStep.c

// Script from http://biostar.stackexchange.com/questions/13869/is-there-a-tool-for-converting-a-fixed-step-wig-to-a-variable-step-wig
// Original authors: Pierre Lindenbaum & Nick
// License is unclear but the code is publicly sitting on the web.
// Added to the AUGUSTUS repository by Katharina Hoff on March 26th 2012
// Compile with gcc -O3 -Wall -ofix2var fixed2variableStep.c

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc,char** argv)
{
int span=1;
int start=0;
int step=0;
char line[BUFSIZ];
while(fgets(line,BUFSIZ,stdin)!=NULL)
 {
 if(strncmp(line,"track",5)==0)
    {
     printf("%s",line);
     continue;
    }
 if(strncmp(line,"fixedStep",9)==0)
    {
    span=1;
    fputs("variableStep",stdout);
    char* saveptr=NULL;
    char* token=line;
    char* ptr;
    while((ptr=strtok_r(token," \n\t", &saveptr))!=NULL)
        {
        char* key=ptr;
        char* value=strchr(ptr,'=');
        if(value==NULL)
            {
            token=saveptr;
            continue;
            }
        *value=0;
        value++;
        if(strcmp(key,"chrom")==0)
            {
            printf(" %s=%s",key,value);
            }
        else if(strcmp(key,"span")==0)
            {
            span=atoi(value);
            }       
        else if(strcmp(key,"start")==0)
            {
            start=atoi(value);
            }
        else if(strcmp(key,"step")==0)
            {
            step=atoi(value);
            }
        token=saveptr;
        }
    printf(" span=%d\n",span);
    continue;
    }
 if(strncmp(line,"0",1)!=0)
   {  
   printf("%d\t%s",start,line);
   }
 start+=step;
 }
return 0;
}
