#!/usr/bin/env nextflow

/*
 * Author: Sidra S
 */ 

/*
 * Pipeline parameters
 */
params.input_file = "data/greetings.csv" - csv file: Hello, Bonjour, Holà

process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
* Use a text replace utility to convert the greeting to uppercase with a UNIX one-liner
*/

process convertToUpper {
    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """

}

workflow {
    
    // create a channel for inputs and read in the csv file - first read the csv file with splitCsv() and then flatten to turn array element from splitCsv() into individual elements
    greeting_ch = Channel.fromPath(params.input_file).splitCsv().flatten()
                         

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
}
