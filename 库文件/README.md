## Table of content
* [库文件](#库文件)
* [动态库](#动态库)


## 库文件
库，其实就是把多个源文件（.c文件），经过一定的翻译，然后打包————到最后只提供给我们一个文件。
Linux库文件是可执行的公共代码，是为了减少程序开发的重复劳动，分为静态库(.a)的动态库(.so)。
静态库文件后缀为.a，程序编译时，将静态文件中的代码复制，拷贝到程序生成的可执行文件中。
静态库的优点就是执行程序时不需要外部的函数库的支持，缺点是如果静态函数库变了，程序必须重新编译。
动态库文件后缀为.so (Shared Object)，动态库在编译是没有被编译进目标代码，而是程序执行到相关代码时才调用库中的函数。


### 动态库 
- .so 结尾
```bash
$ # 查看某个可执行文件的动态库依赖
$ ldd /bin/ls
linux-vdso.so.1 (0x00007fff0e4d5000)
libselinux.so.1 => /usr/lib64/libselinux.so.1 (0x00007fcaf947e000)
libcap.so.2 => /usr/lib64/libcap.so.2 (0x00007fcaf9474000)
libc.so.6 => /usr/lib64/libc.so.6 (0x00007fcaf926b000)
libpcre2-8.so.0 => /usr/lib64/libpcre2-8.so.0 (0x00007fcaf91cf000)
/lib64/ld-linux-x86-64.so.2 (0x00007fcaf94d1000)
```
- so文件的命名规则
```bash
$ libname.so.x.y.z
$  \___/     | | |
$    |       | | release version
$  libname   | minor version
$      major version
$ # major version表示重大升级，不同major version之间的库是不兼容的。
$ # minor version表示增量更新，一般是增加了一些新接口，原来的接口不变。
$ # release version表示库的一些bug修复，性能改进等，不添加任何新的接口，不改变原来的接口。
```

- 动态库查找
```bash
$ # 动态链接库的查找先后顺序为：
$ # (1) LD_LIBRARY_PATH环境变量中的路径
$ # (2) /etc/ld.so.cache缓存文件
$ # (3) /usr/lib和/lib
```

- gcc 是GNU Compiler Collection，原名为Gun C语言编译器，因为它原本只能处理C语言，但gcc很快地扩展，包含很多编译器（C、C++、Objective-C、Ada、Fortran、 Java），可以说gcc是GNU编译器集合；
- g++既可以处理C/C++语言，而gcc只能处理C语言；一般我们使用g++即可；
- gcc/g++就是将包含了代码的文本文件编译（预处理、编译、汇编、链接）成可执行的文件。
- gcc编译选项
使用gcc编译链接时，有两个参数需要注意，一个是-l（小写的L），一个是-L（大写的L）。
Linux有个约定速成的规则，假如库名是name，那么动态链接库文件名就是libname.so。在使用GCC编译链接时，-lname来告诉GCC使用哪个库。
链接时，GCC的链接器ld就会前往LD_LIBRARY_PATH环境变量、/etc/ld.so.cache缓存文件和/usr/lib和/lib目录下去查找libname.so。
我们也可以用-L/path/to/library的方式，让链接器ld去/path/to/library路径下去找库文件。
```bash
gcc -L/path/to/library -lname myfile.c
```

- zlib 安装
```bash
$ sudo yum install zlib-devel -y
```
