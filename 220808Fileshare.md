# How to setup a local-network file-share
The first thing is to evaluate which filesystem structure will hold
our data.
We have the following options:

+ NTFS (Windows)
+ exFAT (Linux)

For both cases we will use a Samba server to share our files. Samba is
the implementation of the SMB protocol.

## Windows Setup
The implementation for the SMB protocol in Windows came with the native folder share
utility that the OS offers. So, Samba Server binaries for Windows it seems that does 
not exist.
One of our needs is to use FOSS (Free and Open Source Software) and it seems that 
the right way to implent a fileshare is to use GNU/Linux. The main drawbacks of 
using Linux as a fileshare are:

- The need to use or exFAT filesystem, that in my opinion is not better that NTFS.
- If you want to use NTFS in Linux, the recommended driver to use the storage is
  ntfs-3g, and it seems not to be completely opensource. And puts some overhead 
  that can corrupt your data.

The reason of I avoided to use a native Windows setup is that your share will
depend on Microsoft's technology stack. Additionally, in Windows it is more
harder to write scripts and automation utilities, since Windows is not a server
oriented operating system. For that reason you will not get natively, rsync, ssh
access or checksuming tools, to name a few.

The only reason to use Windows is its 'reliability', but most of the cloud relies 
on GNU/Linux that for sure will be backed by exFAT filesystems. Anyway, Linux is 
by far another good option. This, at end, does not exclude to you to perform backups
regularly.

Finally in order to have a fileshare easily in you local-network is to use an SFTP
server. In the first stages of this local fileshare, I deployed and the configure
the firewall for a Filezilla SFTP in one of the computers.

## Linux Setup


## References
- [How to share files with Samba](https://www.redhat.com/sysadmin/samba-file-sharing)
