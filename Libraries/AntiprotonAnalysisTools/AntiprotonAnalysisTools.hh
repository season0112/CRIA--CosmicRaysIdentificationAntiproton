#ifndef TimeTools_hh
#define TimeTools_hh

#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <regex>


void Getfilepath(const char *path, const char *filename,  char *filepath)
{
    strcpy(filepath, path); // copy path to filepath (from first until NULL.)
    if(filepath[strlen(path) - 1] != '/') // strlen:calculate the length of string.
        strcat(filepath, "/"); //strcat: add second arguemnt to the end of first argument.
    strcat(filepath, filename); // get the full path of the file.
    //printf("filepath is = %s\n",filepath);
}


bool DeleteFile(const char* path)
{
    printf("\n DeleteFile function start: \n");
    DIR *dir;
    struct dirent *dirinfo;
    struct stat statbuf;  //struct stat: get file status
    char filepath[256] = {0};
    lstat(path, &statbuf); // lstat: get path status to statbuf. return 0 for succee -1 for fail.
    //printf("\n\n\n path now = %s\n",path);
    if (S_ISREG(statbuf.st_mode)) // S_ISREG: interpret the values in a regular file
    {
        remove(path);
        printf("remove successful for %s\n",path);
    }
    else if (S_ISDIR(statbuf.st_mode)) // S_ISDIR: if is a regular directory
    {
        if ((dir = opendir(path)) == NULL){ //opendir:returns a pointer to the directory stream
            printf("Attention: Can not open this dir, please check!\n");
            return 1;}
        while ((dirinfo = readdir(dir)) != NULL) //readdir: read a directory, return a struct dirent pointer.
        {
            Getfilepath(path, dirinfo->d_name, filepath); // 1.Now we have path,dirinfo; this function is to fill filepath.  2.dirinfo is a object of a struct, d_name is a member of this stuct.
            printf("Now the working filepath is: %s\n",filepath);
            if (strcmp(dirinfo->d_name, ".") == 0 || strcmp(dirinfo->d_name, "..") == 0) {// 1.strcmp: compare the two string,  calculate the substract of ASCII of the two first characters, if same then second. 2. ".." is parent  directory, "." is current directory.
                continue;}
            //remove(filepath);
            DeleteFile(filepath);
            rmdir(filepath); //rmdir function:to remove a empty dir.(Attention: rmdir only works if the dir is EMPTY!!!)
            //rmdir(path);
        }
        if ((dirinfo = readdir(dir)) == NULL){
            printf("OK Now the dir is empty!\n");}
        closedir(dir);
    }
    printf("DeleteFile function call ended.\n");
    return 0;
    printf("something else?\n");
}


std::string doubleToString(double price) {
    auto res = std::to_string(price);
    const std::string format("$1");
    try {
        std::regex r("(\\d*)\\.0{6}|");
        std::regex r2("(\\d*\\.{1}0*[^0]+)0*");
        res = std::regex_replace(res, r2, format);
        res = std::regex_replace(res, r, format);
    }
    catch (const std::exception & e) {
        return res;
    }
    return res;
}




#endif 
