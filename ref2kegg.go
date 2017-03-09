package main

import (
	"bufio"
	"log"
	"os"
)

func main() {
	if file, err := os.Open("/export2/home/uesu/db/nr/nr"); err == nil {
		defer file.Close()

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			log.Println(scanner.Text())
		}

		if err = scanner.Err(); err != nil {
			log.Fatal(err)
		}

	} else {
		log.Fatal(err)
	}

}
