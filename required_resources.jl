
# Estimates the resources needed to do a IQA calculation of certain system.
# Assuming that each floating point number occupies 8 Bytes in memory/disk.
# WARNING: do not trust it. It is just a rough approximation. Overflow? BigInt

function required_disk(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp)
  (norbs * (norbs+1)/2 / 1e9) * ( nradial * (lmax+1)^2 + nradial_β * (lmax_β+1)^2) * nat * ncomp * 8 # GB
end

function required_ram(norbs,lmax,nradial,ncomp,size_wfn,nprocs)
  # The maximum required is for the interaction between outside β regions
  # + loading the WFN in nprocs processors
  (norbs * (norbs+1)/2 / 1e9) * nradial * (lmax+1)^2 * 2 * ncomp * 8 + size_wfn * nprocs # GB
end

function io_write(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp)
  required_disk(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp)
end

function io_read(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp,npairs)
  monoread = required_disk(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp)
  monoread + npairs * required_ram(norbs,lmax,nradial,ncomp,0,0) 
end

#------------------------------------------------------------------------------

# System specific
nat      =    12
npairs   = nat*(nat-1)/2
norbs    =  1000
ncomp    =     2 # complex -> 2; non-complex -> 1
size_wfn =    20 # GB

# Parameters chosen
lmax      =   6
lmax_β    =   4
nradial   = 200
nradial_β = 200
nprocs = 10

disk = required_disk(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp)
println("Requires $disk GB of disk space.")
ram = required_ram(norbs,lmax,nradial,ncomp,size_wfn,nprocs)
println("Requires $ram GB of RAM.")
io_w = io_write(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp)
println("Writes $io_w GB to disk (== required).")
io_r = io_read(nat,norbs,lmax,lmax_β,nradial,nradial_β,ncomp,npairs)
println("Reads $io_r GB from disk.")

