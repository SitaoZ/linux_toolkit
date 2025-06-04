## linux 运维
主要收集和linux系统维护与管理相关的命令

- 查看端口
```bash
netstat -tuln
# -t: 显示TCP端口
# -u: 显示UDP端口
# -l: 仅显示监听状态的端口
# -n: 显示数字形式的地址和端口号，而不是尝试解析主机名和服务名

```

- kill listening port
```bash
lsof -i:8000            # 查看该端口的监听情况
kill $(lsof -t -i:8000) # 杀掉该监听端口
```
- terminate called after throwing an instance of 'std::runtime_error'

`terminate called after throwing an instance of 'std::runtime_error'  what():  locale::facet::_S_create_c_locale name not valid`
```bash
$ locale
$ export LC_ALL=C
$ # 最好把它写在.bashrc文件中，这样省的每次都要输入这个命令
$ # LC_ALL=C 是为了去除所有本地化的设置，让命令能正确执行

```


- dnf

DNF 是新一代的rpm软件包管理器。他首先出现在 Fedora 18 这个发行版中。而最近，它取代了yum，正式成为 Fedora 22 的包管理器。
DNF包管理器克服了YUM包管理器的一些瓶颈，提升了包括用户体验，内存占用，依赖分析，运行速度等多方面的内容。DNF使用 RPM, libsolv 和 hawkey 库进行包管理操作。尽管它没有预装在 CentOS7 和 RHEL中，但你可以在使用 YUM 的同时使用 DNF 。
DNF 的最新稳定发行版版本号是 1.0，发行日期是2015年5月11日。 这一版本的额 DNF 包管理器（包括在他之前的所有版本） 都大部分采用 Python 编写，发行许可为GPL v2.

```bash
$ sudo yum install epel-release # yum install epel-release -y
$ sudo yum install dnf

```
