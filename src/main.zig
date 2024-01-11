const std = @import("std");

// entirely based on the python runstats lib by Grant Jenks
// https://github.com/grantjenks/python-runstats

/// Makes a running cummulaive stats struct.
/// ST: The type of the stat variables. Unless you have very speial
/// needs, this should be a float. The significand size should not
/// be larger than your sample size.
pub fn RunningStats(comptime ST: type) type {
    return struct {
        const This = @This();

        n: u64,
        min: ST,
        max: ST,
        eta: ST,
        tau: ST,
        rho: ST,
        phi: ST,

        pub fn new() This {
            return .{
                .n = 0,
                .min = std.math.inf(ST),
                .max = -std.math.inf(ST),
                .eta = 0.0,
                .phi = 0.0,
                .tau = 0.0,
                .rho = 0.0,
            };
        }

        pub fn update(s: *This, x: ST) void {
            s.n += 1;
            s.min = @min(s.min, x);
            s.max = @max(s.max, x);

            const delta = x - s.eta;
            const delta_n = delta / s.n_float();
            const delta_n2 = delta_n * delta_n;
            const term = delta * delta_n * s.n_adj_float(-1);

            // zig fmt: off
            s.eta += delta_n;
            s.phi += term * delta_n2 * (s.n_float() * s.n_float() - 3 * s.n_float() + 3)
                + 6 * delta_n2 * s.rho
                - 4 * delta_n * s.tau;
            s.tau += term * delta_n * s.n_adj_float(-2)
                - 3 * delta_n * s.rho;
            s.rho += term;
        }

        fn n_float(s: *const This) ST {
            return @floatFromInt(s.n);
        }

        fn n_adj_float(s: *const This, d: i64) ST {
            return @floatFromInt(s.n + d);
        }

        pub fn min(s: *const This) ST {
            return s.min;
        }

        pub fn max(s: *const This) ST {
            return s.max;
        }

        pub fn mean(s: *const This) ST {
            return s.eta;
        }

        pub fn variance(s: *const This) ST {
            return s.variance_mdof(1);
        }

        pub fn variance_mdof(s: *const This, dof: u32) ST {
            return s.rho / s.n_float(-@as(i64, dof));
        }

        pub fn stdev(s: *const This) ST {
            return s.stdev_mdof(1);
        }
        pub fn stdev_mdof(s: *const This, dof: u32) ST {
            return @sqrt(s.variance_mdof(dof));
        }

        pub fn skew(s: *const This) ST {
            return @sqrt(s.n_float()) * s.tau / (s.rho * @sqrt(s.rho));
        }

        pub fn kurtosis(s: *const This) ST {
            return s.n_float() * s.phi / (s.rho * s.rho) - 3.0;
        }

        pub fn add(s: *This, o: *const This) void {

            s.max = @max(s.max, o.max);
            s.min = @min(s.min, o.min);

            var tot_n = s.n + o.n;
            var tot_n2 = tot_n * tot_n;
            var tot_n3 = tot_n2 * tot_n;

            var sn2 = s.n * s.n;
            var on2 = o.n * o.n;
            var son = s.n * o.n;

            var d = o.eta - s.eta;
            var d2 = d * d;
            var d3 = d2 * d;
            var d4 = d3 * d;

            // zig fmt: off
            var tot_eta = (s.n * s.eta + o.n * o.eta) / tot_n;

            var tot_rho = s.rho
                + o.rho
                + d2 * s.n * o.n / tot_n;

            var tot_tau = s.tau
                + o.tau
                + d3 * son * (s.n - o.n) / tot_n2
                + 3.0 * d * (s.n * o.rho - o.n * s.rho) / tot_n;


            var tot_phi = s.phi
                + o.phi
                + d4 * son * (sn2 - son + on2) / tot_n3
                + 6.0 * d2 * (sn2 * o.rho + on2 * s.rho) / tot_n2
                + 4.0 * d * (s.n * o.tau - o.n * s.tau) / tot_n;

            s.n = tot_n;
            s.eta = tot_eta;
            s.rho = tot_rho;
            s.tau = tot_tau;
            s.phi = tot_phi;
        }

        pub fn smul(s: *This, scalar: ST) void {
            s.n *= scalar;
            s.phi *= scalar;
            s.tau *= scalar;
            s.rho *= scalar;

        }
    };
}


// --- --- TESTS --- ---
// yeah some day this will happen

test "make it" {
    const Rs16 = RunningStats(f16);
    _ = Rs16.new();
    const Rs32 = RunningStats(f32);
    _ = Rs32.new();
    const Rs64 = RunningStats(f64);
    _ = Rs64.new();
    const Rs80 = RunningStats(f80);
    _ = Rs80.new();
    const Rs128 = RunningStats(f128);
    _ = Rs128.new();
}
