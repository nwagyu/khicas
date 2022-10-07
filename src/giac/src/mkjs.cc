#include <stdio.h>
#include <string.h>

bool translate(const char * filename,FILE * out,FILE * hout){
  char fileext[strlen(filename)+4];
  strcpy(fileext,filename);
  strcpy(fileext+strlen(filename),".js");
  FILE * in=fopen(fileext,"r");
  if (!in)
    return false;
  fprintf(out,"const char %s[]=\"",filename);
  fprintf(hout,"extern const char %s[];\n",filename);
  while (true){
    int c=fgetc(in);
    if (!c || feof(in))
      break;
    if (c=='\n')
      c=' ';
    if (c=='"')
      fputc('\\',out);
    fputc(c,out);
  }
  fclose(in);
  fprintf(out,"%s\n","\";");
  return true;
}

// translate xcasjs.js and qjscalcjs.js to a c file
int main(int argc,const char ** argv){
  FILE * out=fopen("js.c","w");
  FILE *   hout=fopen("js.h","w");
  if (!out || !hout){
    printf("%s\n","Unable to open js.c/js.h for writing");
    return 1;
  }
  fprintf(out,"%s\n","#include \"js.h\"");
  if (argc>1){
    for (int i=1;i<argc;++i){
      translate(argv[i],out,hout);
    }
  }
  else {
    translate("xcasjs",out,hout);
    translate("qjscalcjs",out,hout);
  }
  fclose(out);
  fclose(hout);
  return 0;
}
