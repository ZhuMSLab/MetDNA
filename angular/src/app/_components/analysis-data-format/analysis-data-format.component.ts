import { Component, OnInit } from '@angular/core';

@Component({
  selector: 'app-analysis-data-format',
  templateUrl: './analysis-data-format.component.html',
  styleUrls: ['./analysis-data-format.component.css']
})
export class AnalysisDataFormatComponent implements OnInit {
  public ms1Url: string = '';
  public sampleInfoUrl: string = '';
  constructor() { }

  ngOnInit() {
    this.ms1Url = '/assets/ms1_format.png';
    this.sampleInfoUrl = '/assets/sample_info_format.png';
  }

}
