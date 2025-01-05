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
