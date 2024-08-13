### Slurm


- 显示全部节点
```bash
$ sinfo 
$ sinfo -N
```

- 查看集群的全部节点数
```bash
$ scontrol show node |  grep "mem=" | grep "cpu=" | cut -d , -f 1 | cut -d = -f 3 | awk '{sum+=$1}END{print sum}'
```


- 查看全部任务
```bash
$ squeue
```
