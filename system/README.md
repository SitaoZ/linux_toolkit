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
