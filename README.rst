sharedconfig
============

This is a configuration database to manage hosts that share a common
set of configuration files, but differ in the following ways:

* different machine architectures (e.g. x86_64 and arm64)
* different cpu types (e.g., Skylake and Zen3)
* different client types (e.g., master vs cpu host)
* different machine properties (e.g., number of processors)

These are stored in a way that is consistent with use with
the link farm program ``stow``.  The top-level links are
made to the appropriate arch, cpu, client, and host directories,
then deployed with ``stow *``.  At deployment, links are created
across the whole filesystem.

Updates are done with ``stow -R *``.  Propagation of changes 
to clients can easily be set up with rsync.

