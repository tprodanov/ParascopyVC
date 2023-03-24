
var ploidy = (ploidy === undefined) ? 2 : ploidy;
print('Using ploidy = ' + ploidy)

function record() {
    var sample = SAMPLES[0];
    if (sample.GQ > 0 && sample.GT.split('/').length != ploidy) {
        sample.GQ0 = sample.GQ
        sample.GQ = 0
    }
}
