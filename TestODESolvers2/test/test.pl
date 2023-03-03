#!/usr/bin/perl

use strict;
use warnings;

my %methods = (
    "constant" => {"order" => 0},
    "Euler" => {"order" => 1},
    "RK2" => {"order" => 2},
    "RK3" => {"order" => 3},
    "SSPRK3" => {"order" => 3},
    "RK4" => {"order" => 4},
    "RKF78" => {"order" => 7},
    "DP87" => {"order" => 8},
    # "IMEX-SSP2-222" => {"order" => 3},
    # "IMEX-SSP2-322" => {"order" => 3},
    # "IMEX-SSP2-332" => {"order" => 3},
    # "IMEX-SSP3-332" => {"order" => 3},
    # "IMEX-SSP3-433" => {"order" => 3},
    );

for my $method (keys %methods) {
    my $lines = <<EOF;
ActiveThorns = "
    CarpetX
    IOUtil
    ODESolvers
    TestODESolvers2
"

Cactus::presync_mode = "mixed-error"

CarpetX::xmin = 0
CarpetX::ymin = 0
CarpetX::zmin = 0
CarpetX::xmax = 1
CarpetX::ymax = 1
CarpetX::zmax = 1

CarpetX::ncells_x = 1
CarpetX::ncells_y = 1
CarpetX::ncells_z = 1

CarpetX::blocking_factor_x = 1
CarpetX::blocking_factor_y = 1
CarpetX::blocking_factor_z = 1

CarpetX::ghost_size = 0

CarpetX::dtfac = 1.0
Cactus::cctk_itlast = 1

ODESolvers::method = "$method"

TestODESolvers2::porder = $methods{$method}->{order}

IO::out_dir = \$parfile
IO::parfile_write = "no"
IO::out_every = 1

CarpetX::out_norm_vars = ""

CarpetX::out_tsv_vars = "
    TestODESolvers2::error
    TestODESolvers2::order
    TestODESolvers2::state
"
EOF
    my $fn = "test-$method.par";
    open (my $fh, ">", $fn) or die "Could not open '$fn': $!";
    print $fh $lines;
    close $fh or die "Could not write to '$fn': $!";
}
