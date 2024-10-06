const std = @import("std");

// entirely based on the python runstats lib by Grant Jenks
// https://github.com/grantjenks/python-runstats

/// Makes a running cummulaive stats struct.
/// Sample: The type of the data points and the interior statistics. This is currently
/// required to be a floating point type, and the number of samples should be able to fit
/// into the significand of the float. There is currently no check on this, but there is
/// a plan to add it to debug and release safe builds.
pub fn RunningStats(comptime Sample_: type) type {
    const ti = @typeInfo(Sample_);
    if (ti != .Float) {
        @compileError("Running Stats requires a floating point type for samples");
    }
    return struct {
        const This = @This();
        pub const Sample = Sample_;

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

        pub fn init() This {
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

        pub fn update(this: *This, x: Sample) void {
            this.n += 1;
            this.min = @min(this.min, x);
            this.max = @max(this.max, x);

            const delta = x - this.eta;
            const delta_n = delta / @as(Sample, @floatFromInt(this.n));
            const delta_n2 = delta_n * delta_n;
            const term = delta * delta_n * @as(Sample, @floatFromInt(this.n - 1));

            this.eta += delta_n;
            this.phi += term * delta_n2 * @as(Sample, @floatFromInt((this.n << 1) - 3 * this.n + 3)) +
                6 * delta_n2 * this.rho - 4 * delta_n * this.tau;
            this.tau += term * delta_n * @as(Sample, @floatFromInt(this.n - 2)) -
                3 * delta_n * this.rho;
            this.rho += term;
        }

        pub fn minimum(this: *const This) Sample {
            return this.min;
        }

        pub fn maximum(this: *const This) Sample {
            return this.max;
        }

        pub fn mean(this: *const This) Sample {
            return this.eta;
        }

        pub fn variance(this: *const This) Sample {
            return this.variance_dof(1);
        }

        pub fn variance_dof(this: *const This, dof: u32) Sample {
            return this.rho / @as(Sample, this.n - dof);
        }

        pub fn stdev(this: *const This) Sample {
            return this.stdev_dof(1);
        }

        pub fn stdev_dof(this: *const This, dof: u32) Sample {
            return @sqrt(this.variance_dof(dof));
        }

        pub fn skew(this: *const This) Sample {
            return @sqrt(@as(Sample, @floatFromInt(this.n))) * this.tau / (this.rho * @sqrt(this.rho));
        }

        pub fn kurtosis(this: *const This) Sample {
            return @as(Sample, @floatFromInt(this.n)) * this.phi / (this.rho * this.rho) - 3.0;
        }

        pub fn add(this: *This, o: *const This) void {
            this.max = @max(this.max, o.max);
            this.min = @min(this.min, o.min);

            const tot_n = this.n + o.n;
            const tot_n2 = tot_n * tot_n;
            const tot_n3 = tot_n2 * tot_n;

            const sn2 = this.n * this.n;
            const on2 = o.n * o.n;
            const son = this.n * o.n;

            const d = o.eta - this.eta;
            const d2 = d * d;
            const d3 = d2 * d;
            const d4 = d3 * d;

            // zig fmt: off
            const tot_eta = (this.n * this.eta + o.n * o.eta) / tot_n;

            const tot_rho = this.rho
                + o.rho
                + d2 * this.n * o.n / tot_n;

            const tot_tau = this.tau
                + o.tau
                + d3 * son * (this.n - o.n) / tot_n2
                + 3.0 * d * (this.n * o.rho - o.n * this.rho) / tot_n;


            const tot_phi = this.phi
                + o.phi
                + d4 * son * (sn2 - son + on2) / tot_n3
                + 6.0 * d2 * (sn2 * o.rho + on2 * this.rho) / tot_n2
                + 4.0 * d * (this.n * o.tau - o.n * this.tau) / tot_n;

            this.n = tot_n;
            this.eta = tot_eta;
            this.rho = tot_rho;
            this.tau = tot_tau;
            this.phi = tot_phi;
        }

        pub fn smul(this: *This, scalar: Sample) void {
            this.n *= scalar;
            this.phi *= scalar;
            this.tau *= scalar;
            this.rho *= scalar;

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

        pub fn update(this: *This, x: Sample) void {
            const alpha = 1.0 - this.decay;
            const diff = x - this.mean;
            const incr = alpha * diff;
            this.variance += alpha * (this.decay * diff * diff - this.variance);
            this.mean += incr;
        }

        pub fn mean(this: *const This) Sample {
            return this.mean;
        }

        pub fn variance(this: *const This) Sample {
            return this.variance;
        }

        pub fn stdev(this: *const This) Sample {
            return @sqrt(this.variance);
        }

        pub fn add(this: *This, o: *const This) void {
            this.mean += o.mean;
            this.variance += o.variance;
        }
    };
}

pub fn Corr(SampleT: type) type {
    return struct {
        const This = @This();
        const Sample = SampleT;

        xstatthis:RunningStats(Sample),
        ystatthis:RunningStats(Sample),
        xy: Sample,
        n: i64,

        pub fn update(this: *This, x: Sample, y: Sample) void {
            this.xy += (this.xstatthis.mean() - x)
                * (this.ystatthis.mean() - y)
                * this.n / (this.n + 1);
            this.xstatthis.update(x);
            this.ystatthis.update(y);
            this.n += 1;
        }

        pub fn slope(this: *const This) Sample {
            return this.slope_dof(1.0);
        }

        pub fn slope_dof(this: *const This, dof: Sample) Sample {
            const xx = this.xstatthis.variance_dof(dof) * (this.count - dof);
            return this.xy / xx;
        }

        pub fn intercept(this: *const This) Sample {
            return this.intercept(1.0);
        }

        pub fn intercept_dof(this: *const This, dof: Sample) Sample {
            return this.ystatthis.mean() - this.slope(dof) * this.xstatthis.mean();
        }

        pub fn corr(this: *const This) Sample {
            return this.corr_dof(1.0);
        }

        pub fn corr_dof(this: *const This, dof: Sample) Sample {
            const term = this.xstatthis.stdev_dof(dof) * this.ystatthis.stdev(dof);
            return this.xy / ((this.n - dof) * term);
        }
    };
}


// --- --- tests --- ---
// yeah some day this will happen

const ts = std.testing;

test "make it" {
    const Rs16 = RunningStats(f16);
    var a = Rs16.init();
    a.update(1.0);
    const Rs32 = RunningStats(f32);
    var b = Rs32.init();
    b.update(1.0);
    const Rs64 = RunningStats(f64);
    var c = Rs64.init();
    c.update(1.0);
    const Rs80 = RunningStats(f80);
    var d = Rs80.init();
    d.update(1.0);
    const Rs128 = RunningStats(f128);
    var e = Rs128.init();
    e.update(1.0);

    // _ = RunningStats(i128); // this should fail, not sure how to make it a test
}

test "rstats" {
    const data = [_]f64 {7.0,8.0,9.0,4.0,5.0,6.0,1.0,2.0,2.0,3.0};
    var rs = RunningStats(f64).init();
    for(data) |d| {
        rs.update(d);
    }

    try ts.expectEqual(@as(f64, 1.0), rs.minimum());
    try ts.expectEqual(@as(f64, 9.0), rs.maximum());
    try ts.expectApproxEqAbs(@as(f64, 4.7), rs.mean(), @as(f64, 0.0001));
    //try ts.expectApproxEqAbs(@as(f64, 0.18449), rs.skew(), @as(f64, 0.0001));
    //try ts.expectApproxEqAbs(@as(f64, -1.314017), rs.kurtosis(), @as(f64, 0.0001));
}
