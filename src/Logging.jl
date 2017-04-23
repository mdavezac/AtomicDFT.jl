""" All logging goes through here """
module Log

using DocStringExtensions
import Lumberjack

atomic_dft_log = Lumberjack.LumberMill()

for name in [:debug, :info, :warn, :error]
    @eval begin
        $name(msg::AbstractString; kwargs...) =
            Lumberjack.$name(atomic_dft_log, msg; kwargs...)
        $name(msg::AbstractString, args::Dict) =
            Lumberjack.$name(atomic_dft_log, msg, args)
        $name(msg::AbstractString...) = Lumberjack.$name(atomic_dft_log, msg...)
    end
end
log(truck::Lumberjack.TimberTruck, l::Dict) = log(atomic_dft_log, truck, l)
log(mode::AbstractString) = log(atomic_dft_log, mode)
log(mode::AbstractString, msg::AbstractString)  = log(atomic_dft_log, mode, msg)
log(mode::Symbol, msg::AbstractString)  = log(atomic_dft_log, string(mode), msg)
log(mode::AbstractString, msg::AbstractString, args::Dict) =
    log(atomic_dft_log, mode, msg)
log(mode::Symbol, msg::AbstractString, args::Dict) =
    log(atomic_dft_log, string(mode), msg)
log(mode::AbstractString, args::Dict) = log(atomic_dft_log, mode, args)
log(mode::Symbol, args::Dict) = log(atomic_dft_log, string(mode), args)

configure(; kwargs...) = Lumberjack.configure(atomic_dft_log; kwargs...)

"""
    $(SIGNATURES)

Modifies log-level of all "trucks" in AtomicDFT logs. The input should be one of "debug",
"info", "warn", "error", from least to most verbose.
"""
function set_log_level(level::AbstractString="error")
    for (name, truck) in atomic_dft_log.timber_trucks
        Lumberjack.configure(truck, mode=level)
    end
end
set_log_level(level::Symbol) = set_log_level(string(level))

# Sets log level to error by default
set_log_level()

end
