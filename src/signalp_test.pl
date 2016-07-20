use strict;
use warnings;
use myUtils::SignalPService;

my %secuencias = ('ins' => 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN');
my $signalp_service = myUtils::SignalPService->instance();

$signalp_service->get_signal_peptide_info(\%secuencias);
