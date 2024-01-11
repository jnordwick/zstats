const std = @import("std");

// entirely based on the python runstats lib and the c++ lib
// will give links next checkin

const RunningStats = struct {
    n: u64,
    min: f64,
    max: f64,
    eta: f64,
    tau: f64,
    rho: f64,
    phi: f64,

    pub fn new() RunningStats {
        return RunningStats{
            .n = 0,
            .min = std.math.inf(f64),
            .max = -std.math.inf(f64),
            .eta = 0.0,
            .phi = 0.0,
            .tau = 0.0,
            .rho = 0.0,
        };
    }

    pub fn update(s: *@This(), x: f64) void {
        s.n += 1;
        s.min = @min(s.min, x);
        s.max = @max(s.max, x);

        const delta = x - s.eta;
        const delta_n = delta / s.n_float();
        const delta_n2 = delta_n * delta_n;
        const term = delta * delta_n * s.n_adj_float(-1);

        s.eta += delta_n;

        s.phi += term * delta_n2 * (s.n_float() * s.n_float() - 3 * s.n_float() + 3) + 6 * delta_n2 * s.rho - 4 * delta_n * s.tau;

        s.tau += term * delta_n * s.n_adj_float(-2) - 3 * delta_n * s.rho;

        s.rho += term;
    }

    fn n_float(s: *const @This()) f64 {
        return @floatFromInt(s.n);
    }

    fn n_adj_float(s: *const @This(), d: i64) f64 {
        return @floatFromInt(s.n + d);
    }

    pub fn min(s: *const @This()) f64 {
        return s.min;
    }

    pub fn max(s: *const @This()) f64 {
        return s.max;
    }

    pub fn mean(s: *const @This()) f64 {
        return s.eta;
    }

    pub fn variance(s: *const @This()) f64 {
        return s.variance_mdof(1);
    }

    pub fn variance_mdof(s: *const @This(), dof: u32) f64 {
        return s.rho / s.n_float(-@as(i64, dof));
    }

    pub fn stdev(s: *const @This()) f64 {
        return s.stdev_mdof(1);
    }
    pub fn stdev_mdof(s: *const @This(), dof: u32) f64 {
        return @sqrt(s.variance_mdof(dof));
    }

    pub fn skew(s: *const @This()) f64 {
        return @sqrt(s.n_float()) * s.tau / (s.rho * @sqrt(s.rho));
    }

    pub fn kurtosis(s: *const @This()) f64 {
        return s.n_float() * s.phi / (s.rho * s.rho) - 3.0;
    }
};
