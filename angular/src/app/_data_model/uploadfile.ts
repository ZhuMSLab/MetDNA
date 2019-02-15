export class UploadFile {
    constructor(
        public ms1: File,  // only one file, .csv
        public sampleInfo: File,  // only one file, .csv
        public ms2: File[]  // multiple files, .mgf
      ) { }
}
