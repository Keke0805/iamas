#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#define SIZE_N 1024

void is_dir_exist(const char* path, int* exist){
}
void mkdir_c(const char* path, int* err){
    printf("Starting mkdir = %s\n", path);
    *err = mkdir(path, S_IRWXU | S_IRWXG | S_IRWXO);
    if(err != 0){
        if(errno == EEXIST){
            if(access(path, W_OK) != 0){
                *err = -1;
                printf("no access!\n");
            }else{
                printf("dir exist!\n");
            }
        }
    }else{
        printf("mkdir success!\n");
    }
}
void file_info(const char*filename, int*mode, int *exist, int*time){
    int k;
    int err;
    struct stat buf;
    k=stat(filename,&buf);
    printf("filename = %s\n", filename);
    if(k != 0) {
        *mode = 0; 
        *exist = 0; 
        *time = 0;
        err = mkdir("docs11", S_IRWXU | S_IRWXG | S_IRWXO);
    }else{
        *mode = buf.st_mode;
        if(*mode == 0) *exist=0; else *exist=1;
        *time=buf.st_mtime;
    }
}

void read_array_c_(char *filename, int *arr_len, int *id_arr){
    int flag;
    int total_size;
    char* read_buffer;
    struct stat sstat;
    int in;
    printf("Openfile %s\n", filename);
    stat(filename, &sstat);
    printf("File size: %lld\n", (long long) sstat.st_size);
    read_buffer = (char*)malloc(sstat.st_size+6);
    in = open(filename, O_RDONLY);
    total_size = 0;
    while( (flag = read(in, read_buffer + total_size, SIZE_N)) > 0){
        total_size += flag;
    }
    read_buffer[sstat.st_size+1] = '\0';

    char ch;
    int value = 0;
    int start = 0;
    int j = 0;
    int ch_index;

    char* buffer_index;

    for (buffer_index = read_buffer; buffer_index - read_buffer < total_size; ++buffer_index) {
        if (*buffer_index == '\n') {
            id_arr[j++] = value;
            buffer_index++;
            value = ((*buffer_index) ^ 48);
        } else {
            value = (value << 3) + (value << 1) + ((*buffer_index) ^ 48);
        }
    }

    free(read_buffer);
    close(in);
}
