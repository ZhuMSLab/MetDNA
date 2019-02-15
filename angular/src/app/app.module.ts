import { BrowserModule } from '@angular/platform-browser';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { NgModule } from '@angular/core';
import { HttpModule } from '@angular/http';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';
import { NgUploaderModule } from 'ngx-uploader';
import { ValidateProjectNameDirective } from './validate-project-name.directive';

import { AppComponent } from './app.component';
import { IntroductionComponent, AnalysisComponent,
         HelpComponent, FaqComponent, LogInComponent,
         FileUploadComponent, AnalysisParaSettingComponent,
         AnalysisDataFormatComponent, DemoComponent  } from './_components/index';
import { AppRoutingModule } from './_modules/index';
import { ProjectService, UserService, 
         AuthenticationService, FileService, ParameterService } from './_service/index';


@NgModule({
  declarations: [
    AppComponent,
    IntroductionComponent,
    AnalysisComponent,
    HelpComponent,
    FaqComponent,
    LogInComponent,
    ValidateProjectNameDirective,
    FileUploadComponent,
    AnalysisParaSettingComponent,
    AnalysisDataFormatComponent,
    DemoComponent,
    // FileUploadModule
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    FormsModule,
    ReactiveFormsModule,
    HttpModule,
    NgUploaderModule,
    BrowserAnimationsModule,
  ],
  providers: [ProjectService, UserService, 
              AuthenticationService, FileService,
              ParameterService],  // all services are registered here
  bootstrap: [AppComponent]
})
export class AppModule { }
