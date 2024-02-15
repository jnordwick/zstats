const std = @import("std");

// entirely based on the python runstats lib by Grant Jenks
// https://github.com/grantjenks/python-runstats

/// Makes a running cummulaive stats struct.
/// Sample: The type of the stat variables. Unless you have very speial
/// needs, this should be a float. The significand size should not
/// be larger than your sample size.
pub fn RunningStats(comptime SampleT: type) type {
    return struct {
        const This = @This();
        pub const Sample = SampleT;

        // is mixed with calculations against negatve numbers so just easier tha i32 or casting a lo9t
        /// number of samples
        n: i64,
        /// lowest sample. uses @min for comparisons
        min: Sample,
        /// highest sample. uses @max for comparisons
        max: Sample,
        // these are used to keep running values and from these the moments can be calculatesd
        eta: Sample,
        tau: Sample,
        rho: Sample,
        phi: Sample,

        pub fn new() This {
            return .{
                .n = 0,
                .min = std.math.inf(Sample),
                .max = -std.math.inf(Sample),
                .eta = 0.0,
                .phi = 0.0,
                .tau = 0.0,
                .rho = 0.0,
            };
        }

        pub fn update(s: *This, x: Sample) void {
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

        fn n_float(s: *const This) Sample {
            return @floatFromInt(s.n);
        }

        fn n_adj_float(s: *const This, d: i64) Sample {
            return @floatFromInt(s.n + d);
        }

        pub fn min(s: *const This) Sample {
            return s.min;
        }

        pub fn max(s: *const This) Sample {
            return s.max;
        }

        pub fn mean(s: *const This) Sample {
            return s.eta;
        }

        pub fn variance(s: *const This) Sample {
            return s.variance_mdof(1);
        }

        pub fn variance_mdof(s: *const This, dof: u32) Sample {
            const fdof: f64 = dof;
            return s.rho / (s.n_float() - fdof);
        }

        pub fn stdev(s: *const This) Sample {
            return s.stdev_mdof(1);
        }

        pub fn stdev_mdof(s: *const This, dof: u32) Sample {
            return @sqrt(s.variance_mdof(dof));
        }

        pub fn skew(s: *const This) Sample {
            return @sqrt(s.n_float()) * s.tau / (s.rho * @sqrt(s.rho));
        }

        pub fn kurtosis(s: *const This) Sample {
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

        pub fn smul(s: *This, scalar: Sample) void {
            s.n *= scalar;
            s.phi *= scalar;
            s.tau *= scalar;
            s.rho *= scalar;

        }
    };
}

pub fn ExpStats(SampleT: type) type {
    return struct {
        const This = @This();
        pub const Sample = SampleT;

        /// must be 0 < decay < 1
        decay: Sample = 0.9,
        mean: Sample,
        variance: Sample,

        pub fn init(decay: f64) This {
            return .{.decay = decay, .mean = 0, .variance = 0};
        }

        pub fn update(s: *This, x: Sample) void {
            const alpha = 1.0 - s.decay;
            const diff = x - s.mean;
            const incr = alpha * diff;
            s.variance += alpha * (s.decay * diff * diff - s.variance);
            s.mean += incr;
        }

        pub fn mean(s: *const This) Sample {
            return s.mean;

        }

        pub fn variance(s: *const This) Sample {
            return s.variance;
        }

        pub fn stdev(s: *const This) Sample {
            return @sqrt(s.variance);
        }

        pub fn add(s: *This, o: *const This) void {
            s.mean += o.mean;
            s.variance += o.variance;
        }
    };
}

pub fn Corr(SampleT: type) type {
    return struct {
        const This = @This();
        const Sample = SampleT;

        xstats:RunningStats(Sample),
        ystats:RunningStats(Sample),
        xy: Sample,
        n: i64,

        pub fn update(s: *This, x: Sample, y: Sample) void {
            s.xy += (s.xstats.mean() - x)
                * (s.ystats.mean() - y)
                * s.n / (s.n + 1);
            s.xstats.update(x);
            s.ystats.update(y);
            s.n += 1;
        }

        pub fn slope(s: *const This) Sample {
            return slope_dof(s, 1.0);
        }

        pub fn slope_dof(s: *const This, dof: Sample) Sample {
            const xx = s.xstats.variance_mdof(dof) * (s.count - dof);
            return s.xy / xx;
        }

        pub fn intercept(s: *const This) Sample {
            return intercept(s, 1.0);
        }

        pub fn intercept_dof(s: *const This, dof: Sample) Sample {
            return s.ystats.mean() - s.slope(dof) * s.xstats.mean();
        }

        pub fn corr(s: *const This) Sample {
            return corr_dof(s, 1.0);
        }

        pub fn corr_dof(s: *const This, dof: Sample) Sample {
            const term = s.xstats.stdev_mdof(dof) * s.ystats.stdev(dof);
            return s.xy / ((s.n - dof) * term);
        }
    };
}


// --- --- tests --- ---
// yeah some day this will happen

test "make it" {
    const Rs16 = RunningStats(f16);
    var a = Rs16.new();
    a.update(1.0);
    const Rs32 = RunningStats(f32);
    var b = Rs32.new();
    b.update(1.0);
    const Rs64 = RunningStats(f64);
    var c = Rs64.new();
    c.update(1.0);
    const Rs80 = RunningStats(f80);
    var d = Rs80.new();
    d.update(1.0);
    const Rs128 = RunningStats(f128);
    var e = Rs128.new();
    e.update(1.0);
}
