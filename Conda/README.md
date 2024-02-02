## Table of content
* [创建环境(create)](#创建环境(create))
* 
## Conda tips
Conda软件安装十分便利，可以建立不同的环境对软件进行依赖匹配。

## 创建环境(create)
```bash
$ conda --version                       # 显示conda版本
$ conda create -n env_name python=3.7.2 # 创建环境
$ conda activate env_name               # 激活环境
$ conda deactivate                      # 退出环境
$ conda remove -n env_name --all         # 删除环境

$ conda create --clone env_name1 -n env_name2 #克隆环境到新环境
```

2. 环境(environment) 

```bash
$ conda env list                             # 显示当前所有环境
$ conda info --env                           # 显示当前所有环境
$ conda install -n env_name pkg1 pkg2        # 在指定的环境中安装多个包
$ conda uninstall pkg1 -n env_name           # 在指定的环境下删除包
$ conda env remove -n env_name               # 删除环境
$ conda rename -n env_name env_new_name      # 更改环境名称
$ conda remove -n env_name -c channel_name pkgname # 删除指定环境下，特定来源的包

$ conda env export -n env_name > ENV.yml               # 导出环境信息
$ conda env create -n env_name --file ENV.yml          # 安装指定的环境
```

3. 列出环境中的安装包(list)
```bash
$ conda list                                           # 显示当前环境中已经安装的包
$ conda list --show-channel-urls                       # 列出安装包和来源信息
$ conda list -n env_name                               # 列出env_name环境中的安装包
$ conda list -n env_name --show-channel-urls           # 列出指定环境下安装的全部包和其来源

$ conda list --export > package-list.txt               # 保存安装包便于后续使用
$ conda create -n new_env_name --file package-list.txt # 参考文件重新安装并创建新环境


```

4. 清理安装包的缓存(clean)
```bash
$ conda clean -h
$ conda clean -a # 快速删除
$ conda clean -p # 从可写包缓存中删除没有使用的包，但是不会检查已经安装的包是否有软连接到其中的缓存包
$ conda clean -t # 一键删除anaconda pkgs下面的压缩包

```
5. 包的安装(install)
```bash
$ conda install scipy                         # 当前环境下安装软件
$ conda install -n env_name scipy             # 指定环境下安装软件
$ conda install -c channel_name scipy         # 在指定的环境下安装
$ conda install "pkgname>2.7,<3.5"            # 安装指定的版本
$ conda install "pkgname [version='2.5|3.2']" # 安装指定的版本
$ conda uninstall pkgname
$
$ conda search pkg --info                     # 搜索包的信息
```
6. 配置文件(config)
```bash
$ conda config --show         # 显示已经设置好的配置文件的值
$ conda config --describe     # 显示所有可用的配置文件选项
$ conda config --get channels # 显示已经添加的channel
$ conda config add channels x # 添加镜像
$ # 国内提供conda镜像的大学
$ # 清华大学: https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
$ # 北京外国语大学: https://mirrors.bfsu.edu.cn/help/anaconda/
$ # 南京大学: http://mirrors.nju.edu.cn/
$ # 上海交通大学: https://mirror.sjtu.edu.cn/
$ # 哈尔滨工业大学: http://mirrors.hit.edu.cn/#/home
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
$ vi ~/.condarc # 直接添加镜像网址也可以
```
7.问题
```bash
conda activate bedtools 

CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
To initialize your shell, run

    $ conda init <SHELL_NAME>

Currently supported shells are:
  - bash
  - fish
  - tcsh
  - xonsh
  - zsh
  - powershell

See 'conda init --help' for more information and options.

IMPORTANT: You may need to close and restart your shell after running 'conda init'.
```
```bash
source activate   # 激活环境

conda activate bedtools # 再次激活环境
```
