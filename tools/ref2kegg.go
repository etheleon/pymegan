package main

import (
	//"compress/gzip"
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"regexp"
	"time"
	//"runtime"
	"strconv"
	"strings"
	//"sync"
)

func giWPGI(nrpath string) map[int]int {
	m := make(map[int]int)
	nr, err := os.Open(nrpath)
	if err != nil {
		log.Fatal(err)
	}
	defer nr.Close()

	reader := bufio.NewReader(nr)
	r, _ := regexp.Compile("gi\\|\\d+\\|ref\\|[^|]+")
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		isHeader, err := regexp.Match(">", line[:1])
		if err != nil {
			panic(err)
		}
		if isHeader {
			matches := r.FindAllString(string(line), -1)
			switch {
			case len(matches) > 1:
				nrEntry, children := matches[0], matches[1:]

				parent := strings.Split(nrEntry, "|")
				refseqGI, _ := strconv.Atoi(parent[1])
				m[refseqGI] = refseqGI
				for _, element := range children {
					results := strings.Split(element, "|")
					gi, _ := strconv.Atoi(results[1])
					m[gi] = refseqGI
				}
			case len(matches) == 1:
				parent := strings.Split(matches[0], "|")
				refseqGI, _ := strconv.Atoi(parent[1])
				m[refseqGI] = refseqGI
			}
		}
	}
	return m
}

func giKEGG(filepath string) map[int]string {
	//defer wg.Done()
	m := make(map[int]string)
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatal(err)
	}
	giRegex, _ := regexp.Compile("ncbi-gi:(\\d+)")

	defer file.Close()
	reader := bufio.NewReader(file)
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		row := strings.Split(string(line), "\t")
		matches := giRegex.FindStringSubmatch(row[1])
		gi, _ := strconv.Atoi(matches[1])
		m[gi] = row[0]
	}
	return m
}

func KEGGko(filepath string) map[string]int {
	m := make(map[string]int)
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatal(err)
	}
	koRegex, _ := regexp.Compile("ko:K(\\d+)")

	defer file.Close()
	reader := bufio.NewReader(file)
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		row := strings.Split(string(line), "\t")

		matches := koRegex.FindStringSubmatch(row[1])
		ko, _ := strconv.Atoi(matches[1])
		m[row[0]] = ko
	}
	return m
}

func main() {

	start := time.Now()

	flag.Parse()
	nrpath := flag.Arg(0)
	linker_keggGene2gi := flag.Arg(1)
	linker_koKeggGene := flag.Arg(2)

	giWPGI_map := giWPGI(nrpath)
	giKEGG_map := giKEGG(linker_keggGene2gi)
	KEGGko_map := KEGGko(linker_koKeggGene)

	//Link 3 dicts (KEGGko::ko <- KEGG geneID) <-> (kegggi::KEGG geneID <- gi) <-> (giwp::gi <- WP_refseq's gi)

	keggGI_map := make(map[string]int)
	for gi, geneID := range giKEGG_map {
		//fmt.Printf("%s\t%d\n", geneID, giWPGI_map[gi])
		if WPGI, ok := giWPGI_map[gi]; ok {
			keggGI_map[geneID] = WPGI
		}
	}

	wpgiKO_map := make(map[int]int)

	for geneID, ko := range KEGGko_map {
		//fmt.Printf("%s\t%d\n", geneID, keggGI_map[geneID])
		if gi, ok := keggGI_map[geneID]; ok {
			wpgiKO_map[gi] = ko
		}
	}
	for wpgi, ko := range wpgiKO_map {
		fmt.Printf("%d\t%d\n", wpgi, ko)
	}
	elapsed := time.Since(start)
	log.Printf("Binomial took %s", elapsed)

}
